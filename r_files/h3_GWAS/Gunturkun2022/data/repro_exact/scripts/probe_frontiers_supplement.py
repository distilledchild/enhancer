#!/usr/bin/env python3
"""Probe likely Frontiers supplementary-file URLs and save the first hit."""

from __future__ import annotations

import hashlib
import json
import re
import ssl
import sys
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path


ARTICLE_ID = "790566"
DOI = "10.3389/fpsyt.2022.790566"
JOURNAL_ARTICLE = "fpsyt-13-790566"
SUPP_NAME = "Data_Sheet_1.pdf"


def find_urls_in_html(html: str) -> list[str]:
    urls = set()
    for match in re.finditer(r"https?://[^\"'<>\\\s]+", html):
        url = match.group(0).replace("\\u002F", "/")
        if any(key.lower() in url.lower() for key in ("data_sheet", "supp", "790566", "article/file", "file")):
            urls.add(url)
    for match in re.finditer(r"(?:href|src)=[\"']([^\"']+)[\"']", html, re.I):
        value = match.group(1).replace("\\u002F", "/")
        if any(key.lower() in value.lower() for key in ("data_sheet", "supp", "790566", "article/file", "file")):
            urls.add(urllib.parse.urljoin("https://www.frontiersin.org/", value))
    return sorted(urls)


def candidates(html_path: Path | None = None) -> list[str]:
    base_paths = [
        f"https://www.frontiersin.org/articles/{DOI}/file/{SUPP_NAME}",
        f"https://www.frontiersin.org/articles/{DOI}/file/Data_Sheet_1.PDF",
        f"https://www.frontiersin.org/articles/{DOI}/full#supplementary-material",
        f"https://www.frontiersin.org/files/Articles/{ARTICLE_ID}/{JOURNAL_ARTICLE}-Data_Sheet_1.pdf",
        f"https://www.frontiersin.org/files/Articles/{ARTICLE_ID}/{JOURNAL_ARTICLE}-Data_Sheet_1.PDF",
        f"https://www.frontiersin.org/files/Articles/{ARTICLE_ID}/{JOURNAL_ARTICLE}-s001.pdf",
        f"https://www.frontiersin.org/files/Articles/{ARTICLE_ID}/{JOURNAL_ARTICLE}-supplementary-material.pdf",
        f"https://www.frontiersin.org/files/Articles/{ARTICLE_ID}/{JOURNAL_ARTICLE}-supplementary.zip",
        f"https://www.frontiersin.org/files/Articles/{ARTICLE_ID}/{JOURNAL_ARTICLE}-Data_Sheet_1.docx",
        f"https://www.frontiersin.org/files/Articles/{ARTICLE_ID}/{JOURNAL_ARTICLE}-supplementary-material.zip",
        f"https://www.frontiersin.org/api/articles/{ARTICLE_ID}",
        f"https://www.frontiersin.org/api/articles/10.3389/fpsyt.2022.790566",
    ]
    urls = list(base_paths)
    if html_path and html_path.exists():
        html = html_path.read_text(encoding="utf-8", errors="replace")
        urls.extend(find_urls_in_html(html))
    seen = set()
    unique = []
    for url in urls:
        if url not in seen:
            seen.add(url)
            unique.append(url)
    return unique


def fetch(url: str) -> tuple[int | None, dict[str, str], bytes, str | None]:
    context = ssl._create_unverified_context()
    request = urllib.request.Request(url, headers={"User-Agent": "Codex reproducibility audit"})
    try:
        with urllib.request.urlopen(request, timeout=30, context=context) as response:
            headers = {key.lower(): value for key, value in response.headers.items()}
            body = response.read()
            return response.status, headers, body, None
    except urllib.error.HTTPError as error:
        body = error.read()
        headers = {key.lower(): value for key, value in error.headers.items()}
        return error.code, headers, body, None
    except Exception as error:  # noqa: BLE001
        return None, {}, b"", repr(error)


def extension_from_headers(url: str, headers: dict[str, str]) -> str:
    content_type = headers.get("content-type", "").split(";", 1)[0].lower()
    if "pdf" in content_type or url.lower().endswith(".pdf"):
        return ".pdf"
    if "zip" in content_type or url.lower().endswith(".zip"):
        return ".zip"
    if "json" in content_type:
        return ".json"
    if "html" in content_type:
        return ".html"
    return Path(urllib.parse.urlparse(url).path).suffix or ".body"


def main() -> int:
    if len(sys.argv) != 3:
        print("Usage: probe_frontiers_supplement.py FRONTIERS_FULL_HTML OUTPUT_DIR", file=sys.stderr)
        return 2

    html_path = Path(sys.argv[1])
    out_dir = Path(sys.argv[2])
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    hit_count = 0
    for idx, url in enumerate(candidates(html_path), start=1):
        status, headers, body, error = fetch(url)
        content_type = headers.get("content-type", "")
        sha = hashlib.sha256(body).hexdigest() if body else ""
        row = {
            "idx": idx,
            "url": url,
            "status": status,
            "content_type": content_type,
            "bytes": len(body),
            "sha256": sha,
            "error": error or "",
        }
        rows.append(row)
        useful = bool(status and 200 <= status < 300 and body and b"not found" not in body[:500].lower())
        if useful and (
            "application/pdf" in content_type.lower()
            or "application/zip" in content_type.lower()
            or "json" in content_type.lower()
            or SUPP_NAME.lower() in url.lower()
        ):
            hit_count += 1
            ext = extension_from_headers(url, headers)
            safe_name = f"frontiers_supplement_hit_{hit_count}{ext}"
            (out_dir / safe_name).write_bytes(body)
            (out_dir / f"{safe_name}.url.txt").write_text(url + "\n", encoding="utf-8")

    manifest = out_dir / "frontiers_supplement_probe_manifest.tsv"
    fields = ["idx", "status", "content_type", "bytes", "sha256", "url", "error"]
    with manifest.open("w", encoding="utf-8") as handle:
        handle.write("\t".join(fields) + "\n")
        for row in rows:
            handle.write("\t".join(str(row.get(field, "")) for field in fields) + "\n")

    print(f"probed\t{len(rows)}")
    print(f"useful_hits_saved\t{hit_count}")
    print(f"manifest\t{manifest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
