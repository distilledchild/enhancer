#!/usr/bin/env python3
"""Extract reproducibility-critical tables and methods text from Frontiers XML."""

from __future__ import annotations

import csv
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path


def local_name(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def clean_text(text: str) -> str:
    return re.sub(r"\s+", " ", text).strip()


def text_content(element: ET.Element) -> str:
    return clean_text("".join(element.itertext()))


def find_first(root: ET.Element, tag: str, **attrs: str) -> ET.Element | None:
    for element in root.iter():
        if local_name(element.tag) != tag:
            continue
        if all(element.attrib.get(key) == value for key, value in attrs.items()):
            return element
    return None


def table_rows(table_wrap: ET.Element) -> list[list[str]]:
    rows: list[list[str]] = []
    for tr in table_wrap.iter():
        if local_name(tr.tag) != "tr":
            continue
        cells = [
            text_content(cell)
            for cell in list(tr)
            if local_name(cell.tag) in {"th", "td"}
        ]
        if cells:
            rows.append(cells)
    return rows


def write_tsv(rows: list[list[str]], path: Path) -> None:
    max_cols = max((len(row) for row in rows), default=0)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        for row in rows:
            writer.writerow(row + [""] * (max_cols - len(row)))


def extract_methods(root: ET.Element) -> list[str]:
    method_sections: list[str] = []
    for sec in root.iter():
        if local_name(sec.tag) != "sec":
            continue
        title = ""
        for child in list(sec):
            if local_name(child.tag) == "title":
                title = text_content(child)
                break
        if not title:
            continue
        if any(key in title.lower() for key in ("material", "method", "genotyp", "phenotyp", "association")):
            paragraphs = [
                text_content(child)
                for child in list(sec)
                if local_name(child.tag) == "p" and text_content(child)
            ]
            nested = []
            for child in list(sec):
                if local_name(child.tag) != "sec":
                    continue
                nested_title = ""
                for grand in list(child):
                    if local_name(grand.tag) == "title":
                        nested_title = text_content(grand)
                        break
                nested_paragraphs = [
                    text_content(grand)
                    for grand in list(child)
                    if local_name(grand.tag) == "p" and text_content(grand)
                ]
                if nested_title or nested_paragraphs:
                    nested.append("\n".join([f"## {nested_title}"] + nested_paragraphs))
            block = "\n".join([f"# {title}"] + paragraphs + nested)
            method_sections.append(block)
    return method_sections


def extract_supplements(root: ET.Element) -> list[dict[str, str]]:
    supplements: list[dict[str, str]] = []
    for element in root.iter():
        if local_name(element.tag) != "supplementary-material":
            continue
        row = {"id": element.attrib.get("id", "")}
        for key, value in element.attrib.items():
            row[local_name(key)] = value
        label = find_first(element, "label")
        caption = find_first(element, "caption")
        row["label"] = text_content(label) if label is not None else ""
        row["caption"] = text_content(caption) if caption is not None else ""
        supplements.append(row)
    return supplements


def main() -> int:
    if len(sys.argv) != 3:
        print("Usage: extract_frontiers_xml.py ARTICLE.xml OUTPUT_DIR", file=sys.stderr)
        return 2

    xml_path = Path(sys.argv[1])
    out_dir = Path(sys.argv[2])
    out_dir.mkdir(parents=True, exist_ok=True)

    root = ET.parse(xml_path).getroot()

    summary_lines = [f"source_xml\t{xml_path}"]
    for table_id in ("T1", "T2", "T3"):
        table = find_first(root, "table-wrap", id=table_id)
        if table is None:
            summary_lines.append(f"{table_id}\tmissing")
            continue
        rows = table_rows(table)
        title = ""
        title_el = find_first(table, "label")
        caption_el = find_first(table, "caption")
        if title_el is not None:
            title = text_content(title_el)
        if caption_el is not None:
            title = clean_text(f"{title} {text_content(caption_el)}")
        write_tsv(rows, out_dir / f"{table_id}.tsv")
        (out_dir / f"{table_id}.title.txt").write_text(title + "\n", encoding="utf-8")
        summary_lines.append(f"{table_id}\trows={len(rows)}\tcols={max((len(r) for r in rows), default=0)}\t{title}")

    methods = extract_methods(root)
    (out_dir / "methods_key_points.txt").write_text("\n\n".join(methods) + "\n", encoding="utf-8")
    summary_lines.append(f"methods_sections\t{len(methods)}")

    supplements = extract_supplements(root)
    supp_fields = sorted({key for row in supplements for key in row})
    with (out_dir / "supplementary_materials.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=supp_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(supplements)
    summary_lines.append(f"supplementary_materials\t{len(supplements)}")

    (out_dir / "frontiers_xml_extract_summary.txt").write_text("\n".join(summary_lines) + "\n", encoding="utf-8")
    print("\n".join(summary_lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
