#!/usr/bin/env python3
"""Download GeneNetwork HSNIH-Palmer trait metadata and selected sample data."""

from __future__ import annotations

import concurrent.futures as cf
import csv
import json
import re
import ssl
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any


BASE = "https://genenetwork.org/api/v_pre1"
GROUP = "HSNIH-Palmer"
DATASET = "HSNIH-PalmerPublish"


def trait_ids_from_csv(path: Path) -> list[str]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle)
        header = next(reader)
    ids = []
    for column in header:
        match = re.fullmatch(r"HSR_(\d+)", column)
        if match:
            ids.append(match.group(1))
    return ids


def fetch_json(url: str, retries: int = 3) -> Any:
    context = ssl._create_unverified_context()
    last_error: Exception | None = None
    for attempt in range(retries):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": "Codex reproducibility audit"})
            with urllib.request.urlopen(request, timeout=30, context=context) as response:
                return json.load(response)
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as error:
            last_error = error
            time.sleep(0.5 * (attempt + 1))
    raise RuntimeError(f"failed after {retries} retries: {last_error}")


def fetch_metadata(trait_id: str) -> dict[str, Any]:
    url = f"{BASE}/dataset/{GROUP}/{trait_id}"
    try:
        payload = fetch_json(url)
        if isinstance(payload, dict):
            return {"trait_id": trait_id, "url": url, "ok": True, **payload}
        return {"trait_id": trait_id, "url": url, "ok": False, "error": f"unexpected payload: {type(payload).__name__}"}
    except Exception as error:  # noqa: BLE001 - we want an audit trail, not failure.
        return {"trait_id": trait_id, "url": url, "ok": False, "error": repr(error)}


def write_metadata(rows: list[dict[str, Any]], out_dir: Path) -> None:
    jsonl = out_dir / "HSNIH-PalmerPublish_trait_metadata.jsonl"
    with jsonl.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row, ensure_ascii=False, sort_keys=True) + "\n")

    fields = ["trait_id", "ok", "id", "name", "description", "dataset_type", "url", "error"]
    with (out_dir / "HSNIH-PalmerPublish_trait_metadata.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def select_sample_traits(rows: list[dict[str, Any]]) -> list[str]:
    paper_range = {str(i) for i in range(10418, 10447)}
    keywords = (
        "open field",
        "novel object",
        "social interaction",
        "totaldistance",
        "center",
        "object zone",
        "social zone",
        "latency",
        "frequency",
        "duration",
        "distance",
        "age",
        "batch",
        "coat",
        "sex",
        "gender",
        "center lab",
    )
    selected: set[str] = set(paper_range)
    for row in rows:
        haystack = f"{row.get('name', '')} {row.get('description', '')}".lower()
        if any(keyword in haystack for keyword in keywords):
            selected.add(str(row["trait_id"]))
    return sorted(selected, key=int)


def fetch_sample_data(trait_id: str, out_dir: Path) -> dict[str, Any]:
    url = f"{BASE}/sample_data/{DATASET}/{trait_id}"
    out_path = out_dir / f"{trait_id}.json"
    try:
        payload = fetch_json(url)
        out_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        n = len(payload) if isinstance(payload, list) else None
        nonmissing = None
        if isinstance(payload, list):
            nonmissing = sum(1 for row in payload if row.get("value") not in (None, ""))
        return {"trait_id": trait_id, "url": url, "ok": True, "rows": n, "nonmissing": nonmissing, "path": str(out_path)}
    except Exception as error:  # noqa: BLE001
        err_path = out_dir / f"{trait_id}.error.txt"
        err_path.write_text(f"{url}\n{repr(error)}\n", encoding="utf-8")
        return {"trait_id": trait_id, "url": url, "ok": False, "error": repr(error), "path": str(err_path)}


def write_sample_manifest(rows: list[dict[str, Any]], out_dir: Path) -> None:
    fields = ["trait_id", "ok", "rows", "nonmissing", "url", "path", "error"]
    with (out_dir / "sample_data_manifest.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_keyword_views(metadata_rows: list[dict[str, Any]], out_dir: Path) -> None:
    views = {
        "paper_behavior_traits": re.compile(r"open field|novel object|social interaction|object zone|social zone|center zone", re.I),
        "covariate_candidates": re.compile(r"\bage\b|batch|coat|sex|gender|center", re.I),
    }
    fields = ["trait_id", "name", "description"]
    for label, pattern in views.items():
        matches = [
            {"trait_id": row.get("trait_id", ""), "name": row.get("name", ""), "description": row.get("description", "")}
            for row in metadata_rows
            if row.get("ok") and pattern.search(f"{row.get('name', '')} {row.get('description', '')}")
        ]
        with (out_dir / f"{label}.tsv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            writer.writerows(matches)


def main() -> int:
    if len(sys.argv) != 3:
        print("Usage: scrape_genenetwork.py HSNIH-PalmerPublish.csv OUTPUT_DIR", file=sys.stderr)
        return 2

    csv_path = Path(sys.argv[1])
    out_dir = Path(sys.argv[2])
    out_dir.mkdir(parents=True, exist_ok=True)

    trait_ids = trait_ids_from_csv(csv_path)
    with cf.ThreadPoolExecutor(max_workers=8) as pool:
        metadata_rows = list(pool.map(fetch_metadata, trait_ids))
    metadata_rows.sort(key=lambda row: int(row["trait_id"]))
    write_metadata(metadata_rows, out_dir)
    write_keyword_views(metadata_rows, out_dir)

    selected = select_sample_traits(metadata_rows)
    sample_dir = out_dir / "sample_data_traits"
    sample_dir.mkdir(exist_ok=True)
    (out_dir / "selected_sample_trait_ids.txt").write_text("\n".join(selected) + "\n", encoding="utf-8")
    with cf.ThreadPoolExecutor(max_workers=6) as pool:
        sample_rows = list(pool.map(lambda tid: fetch_sample_data(tid, sample_dir), selected))
    sample_rows.sort(key=lambda row: int(row["trait_id"]))
    write_sample_manifest(sample_rows, out_dir)

    summary = [
        f"csv\t{csv_path}",
        f"traits_in_csv\t{len(trait_ids)}",
        f"metadata_ok\t{sum(1 for row in metadata_rows if row.get('ok'))}",
        f"metadata_failed\t{sum(1 for row in metadata_rows if not row.get('ok'))}",
        f"sample_traits_selected\t{len(selected)}",
        f"sample_data_ok\t{sum(1 for row in sample_rows if row.get('ok'))}",
        f"sample_data_failed\t{sum(1 for row in sample_rows if not row.get('ok'))}",
    ]
    (out_dir / "genenetwork_scrape_summary.txt").write_text("\n".join(summary) + "\n", encoding="utf-8")
    print("\n".join(summary))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
