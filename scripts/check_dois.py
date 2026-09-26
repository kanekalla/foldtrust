"""Validate DOI resolution in all markdown files using Crossref API.

Checks that each DOI resolves via Crossref or doi.org redirect (DataCite/Zenodo).
First authors and titles were verified manually (see docs/benchmark/impact.md:157).

Usage:
    python scripts/check_dois.py
"""

import json
import re
import time
import urllib.request
from pathlib import Path
from typing import Dict, List, Optional


def extract_dois_from_file(filepath: Path) -> List[tuple]:
    """Extract all DOI occurrences from markdown file."""
    dois = []

    with open(filepath, "r", encoding="utf-8") as f:
        content = f.read()

    # Extract DOIs (including those in URLs and link text)
    pattern = r"10\.\d{4,9}/(?:[a-zA-Z0-9.\-_;/]|(?:\([a-zA-Z0-9.\-_]+\)))+"

    for match in re.finditer(pattern, content):
        doi = match.group(0)
        doi = doi.rstrip(".,;:!?")

        # Find line number
        line_num = content[: match.start()].count("\n") + 1

        dois.append((doi, filepath, line_num))

    return dois


def query_crossref(doi: str) -> Optional[Dict]:
    """Query Crossref API for DOI metadata."""
    url = f"https://api.crossref.org/works/{doi}"

    try:
        req = urllib.request.Request(
            url,
            headers={"User-Agent": "FoldTrust/1.0 (mailto:kanekalla@users.noreply.github.com)"},
        )

        with urllib.request.urlopen(req, timeout=10) as response:
            data = json.loads(response.read().decode())
            return data.get("message")

    except urllib.error.HTTPError as e:
        if e.code == 404:
            return None
        raise
    except Exception:
        return None


def check_doi_redirect(doi: str) -> bool:
    """Check if DOI redirects successfully (for DataCite/Zenodo)."""
    url = f"https://doi.org/{doi}"

    try:
        req = urllib.request.Request(url, method="HEAD")
        with urllib.request.urlopen(req, timeout=10) as response:
            return response.status in (200, 301, 302, 303, 307, 308)
    except Exception:
        return False


def check_doi_resolution(doi: str) -> Dict:
    """Check if DOI resolves via Crossref or doi.org redirect."""
    result = {
        "doi": doi,
        "valid": False,
        "method": None,
        "error": None,
    }

    # Try Crossref first (covers most DOIs)
    if not doi.startswith("10.5281/"):
        metadata = query_crossref(doi)
        if metadata:
            result["valid"] = True
            result["method"] = "crossref"
            # Store metadata for reporting
            author_list = metadata.get("author", [])
            first_author = author_list[0].get("family", "Unknown") if author_list else "Unknown"
            result["first_author"] = first_author
            result["title"] = metadata.get("title", ["Unknown"])[0]
            return result

    # Fall back to doi.org redirect for DataCite/Zenodo
    if check_doi_redirect(doi):
        result["valid"] = True
        result["method"] = "doi.org-redirect"
        result["first_author"] = "n/a"
        result["title"] = "n/a"
        return result

    result["error"] = "doi-not-found"
    return result


def main():
    print("=" * 70)
    print("DOI Resolution Check")
    print("=" * 70)

    # Find all .md files, excluding data caches
    repo_root = Path(".")
    md_files = []
    for md_file in repo_root.rglob("*.md"):
        if "/_cache/" in str(md_file) or "/data/_cache/" in str(md_file):
            continue
        md_files.append(md_file)

    print(f"\nSearching {len(md_files)} markdown files...")

    # Extract all DOI occurrences
    all_occurrences = []
    for filepath in sorted(md_files):
        dois = extract_dois_from_file(filepath)
        all_occurrences.extend(dois)

    # Group by DOI
    doi_groups = {}
    for doi, filepath, line_num in all_occurrences:
        if doi not in doi_groups:
            doi_groups[doi] = []
        doi_groups[doi].append((filepath, line_num))

    print(f"Found {len(doi_groups)} unique DOIs ({len(all_occurrences)} occurrences)\n")

    # Check each unique DOI once
    all_results = []
    failed_dois = []

    for i, (doi, occurrences) in enumerate(sorted(doi_groups.items()), 1):
        print(f"[{i}/{len(doi_groups)}] Checking {doi} ({len(occurrences)} occurrences)...")

        # Check resolution once per unique DOI
        result = check_doi_resolution(doi)

        # Record result for each occurrence
        for filepath, line_num in occurrences:
            occurrence_result = result.copy()
            occurrence_result["filepath"] = str(filepath)
            occurrence_result["line"] = line_num
            all_results.append(occurrence_result)

        # Show summary for this DOI
        if result["valid"]:
            print(f"  ✓ Resolved via {result['method']}")
        else:
            print(f"  ✗ Failed: {result['error']}")
            failed_dois.append(doi)

        # Rate limit
        if i < len(doi_groups):
            time.sleep(1)

    # Write output
    output_file = Path("benchmarks/outputs/doi_check.txt")
    output_file.parent.mkdir(parents=True, exist_ok=True)

    with open(output_file, "w") as f:
        f.write("DOI Resolution Results\n")
        f.write("=" * 70 + "\n\n")

        for result in all_results:
            f.write(f"DOI: {result['doi']}\n")
            f.write(f"File: {result['filepath']}:{result['line']}\n")
            f.write(f"Method: {result.get('method', 'n/a')}\n")
            if result["valid"]:
                f.write("Status: VALID\n")
                if result.get("first_author"):
                    f.write(f"First author: {result['first_author']}\n")
                if result.get("title"):
                    f.write(f"Title: {result['title']}\n")
            else:
                f.write("Status: FAILED\n")
                f.write(f"Error: {result['error']}\n")
            f.write("\n")

        f.write("=" * 70 + "\n")
        f.write(f"Total unique DOIs: {len(doi_groups)}\n")
        f.write(f"Total occurrences: {len(all_results)}\n")
        f.write(
            f"Valid: {len(all_results) - len(failed_dois) * len(doi_groups[failed_dois[0]]) if failed_dois else len(all_results)}\n"
        )
        f.write(f"Failed: {len(failed_dois)}\n")

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Total unique DOIs: {len(doi_groups)}")
    print(f"Total occurrences: {len(all_results)}")
    print(f"Valid: {sum(1 for r in all_results if r['valid'])}")
    print(f"Failed: {len(failed_dois)}")
    print(f"\nResults written to {output_file}")

    if failed_dois:
        print("\n❌ Failed DOIs:")
        for doi in failed_dois:
            print(f"  - {doi}")
        return 1
    else:
        print("\n✓ All DOIs resolve successfully!")
        return 0


if __name__ == "__main__":
    exit(main())
