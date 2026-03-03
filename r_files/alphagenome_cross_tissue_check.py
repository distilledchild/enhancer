#!/usr/bin/env python3
"""
AlphaGenome Cross-Tissue Checker (27th Meeting Action Items)

This script validates whether AlphaGenome genuinely distinguishes between different
tissues (e.g., Brain vs. Liver vs. Muscle) using the exact same genomic sequence
and the same marker (H3K27ac).

Input: Same A_dual.csv (or any CSV with chr, start0, end, sequence columns).
Output: 
- A stacked plot comparing the H3K27ac tracks for 3 different tissues.
- A printout of the actual tissues found (biosample names) to share with the PI.
"""

import argparse
import csv
import os
import sys
from pathlib import Path

def ensure_alphagenome_runtime():
    script_dir = Path(__file__).resolve().parent
    venv_python = script_dir / ".venv-alphagenome" / "bin" / "python"
    if venv_python.exists() and str(venv_python) != sys.executable:
        os.execv(str(venv_python), [str(venv_python), str(Path(__file__).resolve()), *sys.argv[1:]])
    raise ModuleNotFoundError("Please run this via .venv-alphagenome/bin/python")

def set_csv_field_size_limit():
    size = sys.maxsize
    while size > 0:
        try:
            csv.field_size_limit(size)
            return size
        except OverflowError:
            size //= 10
    raise RuntimeError("Could not set csv.field_size_limit.")

def load_api_key():
    if os.getenv("ALPHAGENOME_API_KEY"):
        return os.getenv("ALPHAGENOME_API_KEY")
    env_path = Path(__file__).resolve().parent / ".env.alphagenome.local"
    if env_path.exists():
        for line in env_path.read_text().splitlines():
            if "=" in line and line.strip().startswith("ALPHAGENOME_API_KEY"):
                return line.split("=", 1)[1].strip().strip('"').strip("'")
    try:
        from alphagenome import colab_utils
        return colab_utils.get_api_key()
    except Exception:
        raise RuntimeError("No API key found. Set ALPHAGENOME_API_KEY.")

def filter_histone_mark_track(chip_histone_tdata, histone_mark):
    if chip_histone_tdata is None: return None
    metadata = getattr(chip_histone_tdata, "metadata", None)
    if metadata is None: return None
    marks = metadata.get("histone_mark")
    if marks is None: return None
    mask = [m == histone_mark for m in marks]
    if not any(mask): return None
    return chip_histone_tdata.filter_tracks(mask)

def main():
    try:
        import matplotlib.pyplot as plt
        from alphagenome.data import genome
        from alphagenome.models import dna_client
        from alphagenome.visualization import plot_components
    except ModuleNotFoundError:
        ensure_alphagenome_runtime()

    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Path to CSV file")
    parser.add_argument("--output-dir", default="alphagenome_cross_tissue_plots", help="Output directory")
    parser.add_argument("--max-regions", type=int, default=3, help="Max regions to test (to save API calls)")
    args = parser.add_argument_group()
    options = parser.parse_args()

    out_dir = Path(options.output_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    # Tissue definitions (Ontology Terms)
    tissues_to_test = [
        {"label": "Forebrain", "terms": ["UBERON:0001890"]},
        {"label": "Liver", "terms": ["UBERON:0002107"]},
        {"label": "Muscle", "terms": ["UBERON:0001630"]},
    ]

    api_key = load_api_key()
    model = dna_client.create(api_key)
    organism = dna_client.Organism.MUS_MUSCULUS

    set_csv_field_size_limit()

    with open(options.input, "r") as f:
        rows = list(csv.DictReader(f))
    
    # We only need one sequence per unique window to check cross-tissue
    unique_windows = {}
    for r in rows:
        win_id = r.get("window_uid") or f"{r.get('chr')}:{r.get('start0')}-{r.get('end')}"
        if win_id not in unique_windows:
            unique_windows[win_id] = r
            if len(unique_windows) >= options.max_regions:
                break

    print(f"Testing {len(unique_windows)} regions for tissue cross-check...")

    found_biosamples_for_pi = set()

    for win_id, row in unique_windows.items():
        print(f"\n--- Testing Region: {win_id} ---")
        seq = row.get("sequence", "").upper()
        if len(seq) < dna_client.SEQUENCE_LENGTH_1MB:
            pad = dna_client.SEQUENCE_LENGTH_1MB - len(seq)
            left = pad // 2
            right = pad - left
            seq = ("N" * left) + seq + ("N" * right)
        elif len(seq) > dna_client.SEQUENCE_LENGTH_1MB:
            trim = len(seq) - dna_client.SEQUENCE_LENGTH_1MB
            left = trim // 2
            seq = seq[left : left + dna_client.SEQUENCE_LENGTH_1MB]

        interval_1mb = genome.Interval(row.get("chr", "chr1"), 0, dna_client.SEQUENCE_LENGTH_1MB)
        
        compiled_components = []

        for tissue in tissues_to_test:
            print(f"Requesting H3K27ac for {tissue['label']} API...")
            try:
                output = model.predict_sequence(
                    sequence=seq,
                    organism=organism,
                    requested_outputs={dna_client.OutputType.CHIP_HISTONE},
                    interval=interval_1mb,
                    ontology_terms=tissue["terms"]
                )
                
                chip_data = getattr(output, "chip_histone", None)
                h3k27ac_track = filter_histone_mark_track(chip_data, "H3K27ac")
                
                if h3k27ac_track and hasattr(h3k27ac_track, "metadata"):
                    meta = h3k27ac_track.metadata
                    if "biosample_name" in meta:
                        biosamp = meta["biosample_name"]
                        if len(biosamp) > 0:
                            found_biosamples_for_pi.update(list(biosamp))
                        
                    compiled_components.append(
                        plot_components.Tracks(
                            tdata=h3k27ac_track,
                            ylabel_template=f"H3K27ac\n{tissue['label']}",
                            filled=True,
                        )
                    )
                else:
                    print(f" -> No H3K27ac track returned for {tissue['label']}")
            except Exception as e:
                print(f" -> API Error for {tissue['label']}: {e}")
        
        # Plot if we got tracks
        if compiled_components:
            plot_path = out_dir / f"cross_tissue_{win_id.replace(':', '_').replace('-', '_')}.png"
            # Plot the central 100kb
            center = dna_client.SEQUENCE_LENGTH_1MB // 2
            zoom_interval = genome.Interval(row.get("chr", "chr1"), center - 50000, center + 50000)
            
            try:
                plot_components.plot(
                    components=compiled_components,
                    interval=zoom_interval,
                    title=f"Cross-Tissue Overlap Check: {win_id}",
                    despine_keep_bottom=True
                )
                plt.savefig(plot_path, dpi=160, bbox_inches="tight")
                plt.close("all")
                print(f"Saved plot: {plot_path}")
            except Exception as e:
                print(f"Failed to plot {win_id}: {e}")
                plt.close("all")

    print("\n=======================================================")
    print("ACTION ITEM COMPLETED: LIST OF SPECIFIC TISSUES FOUND")
    print("Please send the following list to the PI for transparency:")
    print("-------------------------------------------------------")
    for b in sorted(list(found_biosamples_for_pi)):
        print(f" - {b}")
    print("=======================================================\n")

if __name__ == "__main__":
    main()
