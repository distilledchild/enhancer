#!/usr/bin/env python3
"""
Dual-anchor AlphaGenome runner.

Input: A_dual.csv from alphagenome_dual_export.R
Expected columns:
- pair_uid, role (enhancer/promoter)
- chr, start0, end, sequence
- anchor_start0, anchor_end, window_uid

Output layout:
- <output>/enhancer/*.png
- <output>/promoter/*.png
- <output>/combined/*.png   (1 row x 2 columns per pair)
- dual_run_summary.csv
- dual_pair_summary.tsv

Key optimization:
- API call is made once per unique window_uid (or chr:start-end fallback),
  then reused for all rows in that window group (including enhancer/promoter pair rows).
"""

import argparse
import csv
import hashlib
import math
import os
import re
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np


def ensure_alphagenome_runtime():
    script_dir = Path(__file__).resolve().parent
    venv_python = script_dir / ".venv-alphagenome" / "bin" / "python"
    if venv_python.exists() and str(venv_python) != sys.executable:
        os.execv(str(venv_python), [str(venv_python), str(Path(__file__).resolve()), *sys.argv[1:]])

    raise ModuleNotFoundError(
        "No module named 'alphagenome'. Install dependencies with:\n"
        "  /opt/homebrew/bin/python3.12 -m venv r_files/.venv-alphagenome\n"
        "  r_files/.venv-alphagenome/bin/pip install alphagenome matplotlib\n"
    )


def set_csv_field_size_limit():
    size = sys.maxsize
    while size > 0:
        try:
            csv.field_size_limit(size)
            return size
        except OverflowError:
            size //= 10
    raise RuntimeError("Could not set csv.field_size_limit.")


def parse_args():
    parser = argparse.ArgumentParser(description="AlphaGenome dual-anchor checker")
    parser.add_argument("--input", type=str, required=True, help="Path to A_dual.csv")
    parser.add_argument("--output-dir", type=str, default=None, help="Output root dir")
    parser.add_argument("--api-key", type=str, default=None, help="AlphaGenome API key")
    parser.add_argument("--ontology-terms", type=str, default="UBERON:0002037", help="Comma-separated ontology terms")
    parser.add_argument("--organism", type=str, default="MUS_MUSCULUS", help="dna_client.Organism enum")
    parser.add_argument("--zoom-bp", type=int, default=100000, help="Plot zoom window size")
    parser.add_argument("--fallback-box-bp", type=int, default=10000, help="Fallback box size when anchor interval is missing")
    parser.add_argument("--bg-mode", type=str, default="local_flank", choices=["local_flank", "global_outside"], help="Background mode for FC")
    parser.add_argument("--bg-flank-bp", type=int, default=50000, help="Flank size when bg-mode=local_flank")
    parser.add_argument("--h3k27ac-fc-min", type=float, default=1.25, help="H3K27ac threshold")
    parser.add_argument("--h3k4me1-fc-min", type=float, default=1.10, help="H3K4me1 threshold")
    parser.add_argument("--open-fc-min", type=float, default=1.10, help="Open chromatin threshold")
    parser.add_argument("--h3k4me3-fc-max", type=float, default=1.35, help="Promoter-like threshold")
    parser.add_argument("--max-pairs", type=int, default=0, help="Limit number of pairs (0 = all)")
    parser.add_argument("--no-clean-output", action="store_true", help="Do not remove existing png files in role subdirs")
    return parser.parse_args()


def load_api_key(arg_api_key):
    if arg_api_key:
        return arg_api_key

    script_dir = Path(__file__).resolve().parent
    env_candidates = [
        script_dir / ".env.alphagenome.local",
        Path.cwd() / ".env.alphagenome.local",
    ]
    for env_path in env_candidates:
        if env_path.exists():
            for line in env_path.read_text(encoding="utf-8").splitlines():
                line = line.strip()
                if not line or line.startswith("#") or "=" not in line:
                    continue
                k, v = line.split("=", 1)
                if k.strip() == "ALPHAGENOME_API_KEY" and v.strip():
                    os.environ["ALPHAGENOME_API_KEY"] = v.strip().strip('"').strip("'")
                    break
            if os.getenv("ALPHAGENOME_API_KEY"):
                break

    env_key = os.getenv("ALPHAGENOME_API_KEY")
    if env_key:
        return env_key
    try:
        return colab_utils.get_api_key()
    except Exception as exc:
        raise RuntimeError("Could not resolve API key. Set --api-key or ALPHAGENOME_API_KEY.") from exc


def normalize_sequence(seq, target_len):
    seq = (seq or "").upper()
    seq = re.sub(r"\s+", "", seq)
    seq = re.sub(r"[^ACGTN]", "N", seq)
    n = len(seq)
    if n == target_len:
        return seq, "none"
    if n < target_len:
        pad = target_len - n
        left = pad // 2
        right = pad - left
        return ("N" * left) + seq + ("N" * right), f"pad+{pad}"
    trim = n - target_len
    left = trim // 2
    right = trim - left
    end_idx = n - right if right > 0 else n
    return seq[left:end_idx], f"trim-{trim}"


def safe_name(text):
    return re.sub(r"[^A-Za-z0-9._-]+", "_", text)


def short_hash(text):
    return hashlib.md5(str(text).encode("utf-8")).hexdigest()[:12]


def filter_histone_mark_track(chip_histone_tdata, histone_mark):
    if chip_histone_tdata is None:
        return None
    metadata = getattr(chip_histone_tdata, "metadata", None)
    if metadata is None:
        return None
    marks = metadata.get("histone_mark")
    if marks is None:
        return None
    mask = [m == histone_mark for m in marks]
    if not any(mask):
        return None
    return chip_histone_tdata.filter_tracks(mask)


def compute_interval_metrics(track_data, scoring_intervals_bp, bg_mode="global_outside", bg_flank_bp=50000):
    values = getattr(track_data, "values", None)
    if values is None:
        return {
            "n_tracks": 0,
            "resolution": None,
            "region_bins": 0,
            "bg_bins": 0,
            "region_mean": float("nan"),
            "bg_mean": float("nan"),
            "region_fc_over_bg": float("nan"),
            "region_max": float("nan"),
        }

    n_pos = values.shape[0]
    n_tracks = values.shape[1] if len(values.shape) > 1 else 0
    resolution = int(getattr(track_data, "resolution", 1) or 1)
    if n_tracks == 0 or n_pos == 0:
        return {
            "n_tracks": int(n_tracks),
            "resolution": int(resolution),
            "region_bins": 0,
            "bg_bins": 0,
            "region_mean": float("nan"),
            "bg_mean": float("nan"),
            "region_fc_over_bg": float("nan"),
            "region_max": float("nan"),
        }

    mask = np.zeros(n_pos, dtype=bool)
    for s_bp, e_bp in scoring_intervals_bp:
        s_bin = max(0, int(np.floor(s_bp / resolution)))
        e_bin = min(n_pos, int(np.ceil(e_bp / resolution)))
        if e_bin <= s_bin:
            continue
        mask[s_bin:e_bin] = True

    if not mask.any():
        return {
            "n_tracks": int(n_tracks),
            "resolution": int(resolution),
            "region_bins": 0,
            "bg_bins": 0,
            "region_mean": float("nan"),
            "bg_mean": float("nan"),
            "region_fc_over_bg": float("nan"),
            "region_max": float("nan"),
        }

    region_vals = values[mask, :]
    if bg_mode == "local_flank":
        flank_bins = max(1, int(np.ceil(max(1, bg_flank_bp) / resolution)))
        bg_mask = np.zeros(n_pos, dtype=bool)
        for s_bp, e_bp in scoring_intervals_bp:
            s_bin = max(0, int(np.floor(s_bp / resolution)))
            e_bin = min(n_pos, int(np.ceil(e_bp / resolution)))
            if e_bin <= s_bin:
                continue
            left_lo = max(0, s_bin - flank_bins)
            left_hi = s_bin
            right_lo = e_bin
            right_hi = min(n_pos, e_bin + flank_bins)
            if left_hi > left_lo:
                bg_mask[left_lo:left_hi] = True
            if right_hi > right_lo:
                bg_mask[right_lo:right_hi] = True
        bg_mask = bg_mask & (~mask)
        bg_vals = values[bg_mask, :] if bg_mask.any() else values[~mask, :]
    else:
        bg_vals = values[~mask, :]

    region_mean = float(region_vals.mean())
    region_max = float(region_vals.max())
    bg_mean = float(bg_vals.mean()) if bg_vals.size > 0 else float("nan")
    eps = 1e-12
    region_fc_over_bg = float(region_mean / (bg_mean + eps)) if bg_mean == bg_mean else float("nan")
    return {
        "n_tracks": int(n_tracks),
        "resolution": int(resolution),
        "region_bins": int(mask.sum()),
        "bg_bins": int(bg_vals.shape[0]) if bg_vals is not None else 0,
        "region_mean": region_mean,
        "bg_mean": bg_mean,
        "region_fc_over_bg": region_fc_over_bg,
        "region_max": region_max,
    }


def parse_anchor_interval_rel(row, seq_len_required, fallback_bp=10000):
    try:
        ws = int(float(str(row.get("start0", "")).strip()))
        a0 = int(float(str(row.get("anchor_start0", "")).strip()))
        a1 = int(float(str(row.get("anchor_end", "")).strip()))
        rel_s = max(0, a0 - ws)
        rel_e = min(seq_len_required, a1 - ws)
        if rel_e > rel_s:
            return [(rel_s, rel_e)], "anchor_interval"
    except Exception:
        pass

    center = seq_len_required // 2
    half = max(1, fallback_bp // 2)
    return [(max(0, center - half), min(seq_len_required, center + half))], "center_fallback"


def build_anchor_centered_plot_interval(chr_name, seq_len_required, scoring_intervals_rel_bp, zoom_bp):
    if zoom_bp is None or zoom_bp <= 0:
        return genome.Interval(chr_name, 0, seq_len_required)

    if not scoring_intervals_rel_bp:
        center = seq_len_required // 2
    else:
        s_min = min(s for s, _ in scoring_intervals_rel_bp)
        e_max = max(e for _, e in scoring_intervals_rel_bp)
        center = (s_min + e_max) // 2

    half = max(1, zoom_bp // 2)
    start = max(0, center - half)
    end = start + zoom_bp
    if end > seq_len_required:
        end = seq_len_required
        start = max(0, end - zoom_bp)
    if end <= start:
        start, end = 0, seq_len_required
    return genome.Interval(chr_name, int(start), int(end))


def draw_role_boxes(intervals_rel_bp, role):
    fig = plt.gcf()
    if fig is None or not fig.axes:
        return
    if str(role).lower() == "promoter":
        face, edge, alpha = "#f8bbd0", "#e91e63", 0.28
    else:
        face, edge, alpha = "#f1c40f", "#c49b06", 0.24
    for s_bp, e_bp in intervals_rel_bp:
        for ax in fig.axes:
            ax.axvspan(s_bp, e_bp, facecolor=face, edgecolor=edge, linewidth=0.9, alpha=alpha)


def combine_two_plots(enhancer_png, promoter_png, out_png, combined_title):
    img1 = mpimg.imread(enhancer_png)
    img2 = mpimg.imread(promoter_png)
    fig, axes = plt.subplots(1, 2, figsize=(18, 6), dpi=160)
    axes[0].imshow(img1)
    axes[0].axis("off")
    axes[0].set_title("Enhancer anchor (yellow box)")
    axes[1].imshow(img2)
    axes[1].axis("off")
    axes[1].set_title("Promoter anchor (pink box)")
    fig.suptitle(combined_title, fontsize=10)
    plt.tight_layout()
    fig.savefig(out_png, dpi=160, bbox_inches="tight")
    plt.close(fig)


def parse_int(value, default=None):
    try:
        return int(float(str(value)))
    except Exception:
        return default


def bp_to_size_label(bp):
    if bp is None:
        return "NA"
    if bp % 1000 == 0:
        return f"{bp // 1000}k"
    return f"{bp}bp"


def main():
    try:
        import matplotlib.image as _mpimg
        import matplotlib.pyplot as _plt
        from alphagenome import colab_utils as _colab_utils
        from alphagenome.data import genome as _genome
        from alphagenome.models import dna_client as _dna_client
        from alphagenome.visualization import plot_components as _plot_components
    except ModuleNotFoundError as exc:
        if exc.name in {"alphagenome", "matplotlib", "matplotlib.pyplot", "matplotlib.image"}:
            ensure_alphagenome_runtime()
        raise

    global mpimg, plt, colab_utils, genome, dna_client, plot_components
    mpimg = _mpimg
    plt = _plt
    colab_utils = _colab_utils
    genome = _genome
    dna_client = _dna_client
    plot_components = _plot_components

    args = parse_args()
    input_path = Path(args.input).expanduser()
    if not input_path.exists():
        raise FileNotFoundError(f"Input not found: {input_path}")

    output_root = Path(args.output_dir).expanduser() if args.output_dir else (input_path.parent / "alphagenome_dual_plots")
    enhancer_dir = output_root / "enhancer"
    promoter_dir = output_root / "promoter"
    combined_dir = output_root / "combined"
    for d in [output_root, enhancer_dir, promoter_dir, combined_dir]:
        d.mkdir(parents=True, exist_ok=True)

    if not args.no_clean_output:
        for d in [enhancer_dir, promoter_dir, combined_dir]:
            for p in d.glob("*.png"):
                p.unlink(missing_ok=True)

    csv_limit = set_csv_field_size_limit()
    ontology_terms = [x.strip() for x in args.ontology_terms.split(",") if x.strip()]
    api_key = load_api_key(args.api_key)
    model = dna_client.create(api_key)
    try:
        organism = getattr(dna_client.Organism, args.organism)
    except AttributeError as exc:
        raise ValueError(f"Invalid organism: {args.organism}") from exc

    with input_path.open("r", newline="") as f:
        rows = list(csv.DictReader(f))

    required_cols = ["pair_uid", "role", "chr", "start0", "end", "sequence", "anchor_start0", "anchor_end"]
    missing = [c for c in required_cols if c not in (rows[0].keys() if rows else [])]
    if missing:
        raise ValueError(f"Missing required columns in input: {missing}")

    if args.max_pairs > 0:
        keep_pairs = []
        seen = set()
        for r in rows:
            puid = r.get("pair_uid", "")
            if puid not in seen:
                seen.add(puid)
                keep_pairs.append(puid)
            if len(keep_pairs) >= args.max_pairs:
                break
        keep_pairs = set(keep_pairs)
        rows = [r for r in rows if r.get("pair_uid", "") in keep_pairs]

    # Build deterministic filenames based on gene symbol.
    # If one gene has multiple loop pairs, suffix _2, _3... to avoid overwrite.
    pair_meta = {}
    pair_order = []
    for row in rows:
        pair_uid = row.get("pair_uid", "")
        if pair_uid not in pair_meta:
            pair_order.append(pair_uid)
            pair_meta[pair_uid] = {
                "gene_symbol": row.get("gene_symbol", ""),
                "enhancer_anchor_uid": "",
                "promoter_anchor_uid": "",
                "enhancer_size_label": "NA",
                "promoter_size_label": "NA",
            }
        role = (row.get("role") or "").strip().lower()
        a_chr = row.get("anchor_chr", "")
        a0 = parse_int(row.get("anchor_start0"))
        a1 = parse_int(row.get("anchor_end"))
        if a0 is not None and a1 is not None:
            anchor_uid = f"{a_chr}:{a0}-{a1}"
            size_label = bp_to_size_label(max(0, a1 - a0))
        else:
            anchor_uid = row.get("anchor_uid", "")
            size_label = bp_to_size_label(parse_int(row.get("anchor_width_bp")))
        if role == "enhancer":
            pair_meta[pair_uid]["enhancer_anchor_uid"] = anchor_uid
            pair_meta[pair_uid]["enhancer_size_label"] = size_label
        elif role == "promoter":
            pair_meta[pair_uid]["promoter_anchor_uid"] = anchor_uid
            pair_meta[pair_uid]["promoter_size_label"] = size_label

    gene_counts = defaultdict(int)
    pair_stub = {}
    for pair_uid in pair_order:
        gene_symbol = (pair_meta.get(pair_uid, {}).get("gene_symbol") or "NA").strip()
        gene_base = safe_name(gene_symbol) or "NA"
        gene_counts[gene_base] += 1
        idx = gene_counts[gene_base]
        suffix = "" if idx == 1 else f"_{idx}"
        pair_stub[pair_uid] = f"{gene_base}{suffix}"

    seq_len_required = dna_client.SEQUENCE_LENGTH_1MB
    print(f"[Dual] Input: {input_path}")
    print(f"[Dual] Rows: {len(rows)} | csv_field_size_limit={csv_limit}")
    print(f"[Dual] Output: {output_root}")
    print(f"[Dual] Background: mode={args.bg_mode}, flank_bp={args.bg_flank_bp}")

    # Group by window. One API call per unique window.
    groups = defaultdict(list)
    for i, row in enumerate(rows):
        key = row.get("window_uid") or f"{row.get('chr')}:{row.get('start0')}-{row.get('end')}"
        groups[key].append(i)
    print(f"[Dual] Unique windows (API calls target): {len(groups)}")

    summary_rows = []
    pair_role_to_plot = defaultdict(dict)
    t0 = time.time()
    api_calls = 0
    api_errors = 0

    for gi, (window_key, idxs) in enumerate(groups.items(), start=1):
        row0 = rows[idxs[0]]
        chr_name = row0.get("chr") or "chr1"
        raw_seq = row0.get("sequence") or ""
        seq, adjust = normalize_sequence(raw_seq, seq_len_required)
        print(f"[Dual][window {gi}/{len(groups)}] {window_key} | rows_in_window={len(idxs)} | adjust={adjust}")

        try:
            predict_kwargs = {
                "sequence": seq,
                "organism": organism,
                "requested_outputs": {
                    dna_client.OutputType.ATAC,
                    dna_client.OutputType.DNASE,
                    dna_client.OutputType.CHIP_HISTONE,
                },
                "interval": genome.Interval(chr_name, 0, seq_len_required),
                "ontology_terms": ontology_terms,
            }
            output = model.predict_sequence(**predict_kwargs)
            api_calls += 1
        except Exception as exc:
            api_errors += 1
            for idx in idxs:
                row = rows[idx]
                summary_rows.append({
                    "pair_uid": row.get("pair_uid"),
                    "role": row.get("role"),
                    "window_uid": row.get("window_uid") or window_key,
                    "status": "error",
                    "error": str(exc),
                    "plot_png": "",
                })
            print(f"  -> ERROR API window: {exc}")
            continue

        h3k27ac_track = filter_histone_mark_track(output.chip_histone, "H3K27ac")
        h3k4me1_track = filter_histone_mark_track(output.chip_histone, "H3K4me1")
        h3k4me3_track = filter_histone_mark_track(output.chip_histone, "H3K4me3")

        for idx in idxs:
            row = rows[idx]
            role = (row.get("role") or "").strip().lower()
            pair_uid = row.get("pair_uid") or f"pair_{idx+1}"
            loop_id = row.get("loop_id") or ""
            gene_symbol = row.get("gene_symbol") or ""
            anchor_uid = row.get("anchor_uid") or ""
            anchor_chr = row.get("anchor_chr") or chr_name
            anchor_start0 = parse_int(row.get("anchor_start0"))
            anchor_end = parse_int(row.get("anchor_end"))
            if anchor_start0 is not None and anchor_end is not None:
                anchor_uid = f"{anchor_chr}:{anchor_start0}-{anchor_end}"
                anchor_size_label = bp_to_size_label(max(0, anchor_end - anchor_start0))
            else:
                anchor_size_label = bp_to_size_label(parse_int(row.get("anchor_width_bp")))
            scoring_intervals, interval_source = parse_anchor_interval_rel(row, seq_len_required, fallback_bp=args.fallback_box_bp)

            atac_m = compute_interval_metrics(output.atac, scoring_intervals, bg_mode=args.bg_mode, bg_flank_bp=args.bg_flank_bp)
            dnase_m = compute_interval_metrics(output.dnase, scoring_intervals, bg_mode=args.bg_mode, bg_flank_bp=args.bg_flank_bp)
            h3k27_m = compute_interval_metrics(h3k27ac_track, scoring_intervals, bg_mode=args.bg_mode, bg_flank_bp=args.bg_flank_bp)
            h3k4m1_m = compute_interval_metrics(h3k4me1_track, scoring_intervals, bg_mode=args.bg_mode, bg_flank_bp=args.bg_flank_bp)
            h3k4m3_m = compute_interval_metrics(h3k4me3_track, scoring_intervals, bg_mode=args.bg_mode, bg_flank_bp=args.bg_flank_bp)

            if atac_m["n_tracks"] > 0:
                open_src, open_m = "ATAC", atac_m
            elif dnase_m["n_tracks"] > 0:
                open_src, open_m = "DNASE", dnase_m
            else:
                open_src, open_m = "none", compute_interval_metrics(None, scoring_intervals, bg_mode=args.bg_mode, bg_flank_bp=args.bg_flank_bp)

            h3k27_ok = h3k27_m["n_tracks"] > 0 and h3k27_m["region_fc_over_bg"] == h3k27_m["region_fc_over_bg"] and h3k27_m["region_fc_over_bg"] >= args.h3k27ac_fc_min
            h3k4m1_ok = h3k4m1_m["n_tracks"] > 0 and h3k4m1_m["region_fc_over_bg"] == h3k4m1_m["region_fc_over_bg"] and h3k4m1_m["region_fc_over_bg"] >= args.h3k4me1_fc_min
            open_ok = open_m["n_tracks"] > 0 and open_m["region_fc_over_bg"] == open_m["region_fc_over_bg"] and open_m["region_fc_over_bg"] >= args.open_fc_min
            promoter_like = h3k27_ok and h3k4m3_m["n_tracks"] > 0 and h3k4m3_m["region_fc_over_bg"] == h3k4m3_m["region_fc_over_bg"] and h3k4m3_m["region_fc_over_bg"] >= args.h3k4me3_fc_max

            support_count = int(bool(h3k4m1_ok)) + int(bool(open_ok))
            if promoter_like:
                decision_class = "promoter_like"
            elif h3k27_ok and h3k4m1_ok and open_ok:
                decision_class = "high_confidence"
            elif h3k27_ok and (h3k4m1_ok or open_ok):
                decision_class = "enhancer_like"
            else:
                decision_class = "not_enhancer_like"

            # Role-specific plot path.
            role_dir = promoter_dir if role == "promoter" else enhancer_dir
            base = pair_stub.get(pair_uid, safe_name(gene_symbol) or short_hash(pair_uid))
            plot_path = role_dir / f"{base}_{role}.png"

            try:
                components = []
                if open_src == "ATAC":
                    components.append(plot_components.Tracks(tdata=output.atac, ylabel_template="ATAC\n{biosample_name}", filled=True))
                elif open_src == "DNASE":
                    components.append(plot_components.Tracks(tdata=output.dnase, ylabel_template="DNASE\n{biosample_name}", filled=True))
                if h3k27_m["n_tracks"] > 0:
                    components.append(plot_components.Tracks(tdata=h3k27ac_track, ylabel_template="H3K27ac\n{biosample_name}", filled=True))
                if h3k4m1_m["n_tracks"] > 0:
                    components.append(plot_components.Tracks(tdata=h3k4me1_track, ylabel_template="H3K4me1\n{biosample_name}", filled=True))
                if h3k4m3_m["n_tracks"] > 0:
                    components.append(plot_components.Tracks(tdata=h3k4me3_track, ylabel_template="H3K4me3\n{biosample_name}", filled=True))

                if components:
                    zoom_interval = build_anchor_centered_plot_interval(
                        chr_name=chr_name,
                        seq_len_required=seq_len_required,
                        scoring_intervals_rel_bp=scoring_intervals,
                        zoom_bp=args.zoom_bp,
                    )
                    plot_components.plot(
                        components=components,
                        interval=zoom_interval,
                        title=(
                            f"{gene_symbol} | role={role} | anchor={anchor_uid} ({anchor_size_label}) "
                            f"| window={chr_name}:{row.get('start0')}-{row.get('end')}"
                        ),
                        despine_keep_bottom=True,
                    )
                    draw_role_boxes(scoring_intervals, role=role)
                    plt.savefig(plot_path, dpi=160, bbox_inches="tight")
                plt.close("all")
            except Exception:
                plt.close("all")
                plot_path = Path("")

            if str(plot_path):
                pair_role_to_plot[pair_uid][role] = str(plot_path)

            summary_rows.append({
                "pair_uid": pair_uid,
                "gene_symbol": gene_symbol,
                "loop_id": loop_id,
                "role": role,
                "anchor_uid": anchor_uid,
                "window_uid": row.get("window_uid") or window_key,
                "shared_window_used": row.get("shared_window_used", ""),
                "status": "ok",
                "error": "",
                "interval_source": interval_source,
                "bg_mode": args.bg_mode,
                "bg_flank_bp": args.bg_flank_bp,
                "decision_class": decision_class,
                "support_count": support_count,
                "open_source": open_src,
                "h3k27ac_fc": h3k27_m["region_fc_over_bg"],
                "h3k4me1_fc": h3k4m1_m["region_fc_over_bg"],
                "open_fc": open_m["region_fc_over_bg"],
                "h3k4me3_fc": h3k4m3_m["region_fc_over_bg"],
                "plot_png": str(plot_path),
            })

    # Build combined 1x2 images per pair.
    n_combined = 0
    for pair_uid, rp in pair_role_to_plot.items():
        enh_png = rp.get("enhancer")
        prom_png = rp.get("promoter")
        if not enh_png or not prom_png:
            continue
        base = pair_stub.get(pair_uid, safe_name(pair_uid) or short_hash(pair_uid))
        combined_png = combined_dir / f"{base}_combined.png"
        meta = pair_meta.get(pair_uid, {})
        combined_title = (
            f"{meta.get('gene_symbol','NA')} | enhancer={meta.get('enhancer_anchor_uid','NA')} ({meta.get('enhancer_size_label','NA')}) "
            f"| promoter={meta.get('promoter_anchor_uid','NA')} ({meta.get('promoter_size_label','NA')})"
        )
        try:
            combine_two_plots(enh_png, prom_png, combined_png, combined_title)
            n_combined += 1
        except Exception:
            continue

    # Write summary files.
    summary_path = output_root / "dual_run_summary.csv"
    fieldnames = [
        "pair_uid", "gene_symbol", "loop_id", "role", "anchor_uid", "window_uid", "shared_window_used",
        "status", "error", "interval_source", "bg_mode", "bg_flank_bp", "decision_class", "support_count",
        "open_source", "h3k27ac_fc", "h3k4me1_fc", "open_fc", "h3k4me3_fc", "plot_png",
    ]
    with summary_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(summary_rows)

    n_ok = sum(1 for r in summary_rows if r.get("status") == "ok")
    n_err = len(summary_rows) - n_ok
    role_counts = defaultdict(int)
    class_counts = defaultdict(int)
    for r in summary_rows:
        role_counts[r.get("role", "")] += 1
        if r.get("status") == "ok":
            class_counts[r.get("decision_class", "")] += 1

    pair_summary_path = output_root / "dual_pair_summary.tsv"
    with pair_summary_path.open("w", newline="") as f:
        f.write("metric\tvalue\n")
        f.write(f"rows_total\t{len(summary_rows)}\n")
        f.write(f"rows_ok\t{n_ok}\n")
        f.write(f"rows_error\t{n_err}\n")
        f.write(f"unique_pairs\t{len({r.get('pair_uid') for r in summary_rows})}\n")
        f.write(f"role_enhancer_rows\t{role_counts.get('enhancer', 0)}\n")
        f.write(f"role_promoter_rows\t{role_counts.get('promoter', 0)}\n")
        f.write(f"api_calls_unique_windows\t{api_calls}\n")
        f.write(f"api_errors_windows\t{api_errors}\n")
        f.write(f"combined_png_count\t{n_combined}\n")
        for cls in ["high_confidence", "enhancer_like", "promoter_like", "not_enhancer_like"]:
            f.write(f"{cls}_rows\t{class_counts.get(cls, 0)}\n")

    elapsed = time.time() - t0
    print(f"[Dual] Done. rows_ok={n_ok}, rows_error={n_err}, elapsed={elapsed:.1f}s")
    print(f"[Dual] API calls={api_calls} (unique windows), api_errors={api_errors}")
    print(f"[Dual] enhancer_png={len(list(enhancer_dir.glob('*.png')))}, promoter_png={len(list(promoter_dir.glob('*.png')))}, combined_png={len(list(combined_dir.glob('*.png')))}")
    print(f"[Dual] Summary: {summary_path}")
    print(f"[Dual] Pair summary: {pair_summary_path}")


if __name__ == "__main__":
    main()
