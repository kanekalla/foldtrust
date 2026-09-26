"""Validate DOIs in all markdown files using Crossref API.

Checks each DOI by querying the Crossref API and verifying metadata.
Falls back to HEAD request for DataCite/Zenodo DOIs.

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
    """Extract DOIs from markdown file."""
    dois = []

    with open(filepath, "r") as f:
        content = f.read()

    # Extract DOIs - use negative lookahead to stop before markdown closing )
    # DOI suffix can contain alphanumerics, dots, hyphens, underscores, slashes
    # and balanced parentheses (e.g., 10.1016/0888-7543(91)90503-7)
    
    # Basic pattern: 10.xxxx/yyy where yyy doesn't contain ) unless preceded by (
    pattern = r"10\.\d{4,9}/(?:[a-zA-Z0-9.\-_;/]|(?:\([a-zA-Z0-9.\-_]+\)))+"
    
    for match in re.finditer(pattern, content):
        doi = match.group(0)
        # Strip any trailing punctuation
        doi = doi.rstrip(".,;:!?")
        dois.append((doi, filepath))

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
            # Check if we get a 3xx redirect to a non-error page
            return response.status in (200, 301, 302, 303, 307, 308)
    except Exception:
        return False


def format_authors(authors: List[Dict]) -> str:
    """Format author list from Crossref metadata."""
    if not authors:
        return "Unknown"

    first_author = authors[0]
    family = first_author.get("family", "")
    given = first_author.get("given", "")

    if len(authors) == 1:
        return f"{family}, {given}" if given else family
    else:
        return f"{family} et al."


def check_doi(doi: str) -> Dict:
    """Check a single DOI and return validation results."""
    result = {
        "doi": doi,
        "valid": False,
        "error": None,
        "title": None,
        "authors": None,
    }

    # Try Crossref first
    metadata = query_crossref(doi)

    if metadata:
        title = metadata.get("title", [""])[0]
        authors = format_authors(metadata.get("author", []))

        result["valid"] = True
        result["title"] = title
        result["authors"] = authors
        return result

    # Fallback: try HEAD request for DataCite/Zenodo DOIs
    if check_doi_redirect(doi):
        result["valid"] = True
        result["title"] = "(DataCite/Zenodo - redirect OK)"
        result["authors"] = "(not available)"
        return result

    result["error"] = "DOI not found or inaccessible"
    return result


def main():
    print("=" * 70)
    print("DOI Validation")
    print("=" * 70)

    # Find all .md files, excluding data caches
    repo_root = Path(".")
    md_files = []
    for md_file in repo_root.rglob("*.md"):
        # Skip cache directories
        if "/_cache/" in str(md_file) or "/data/_cache/" in str(md_file):
            continue
        md_files.append(md_file)

    print(f"\nSearching {len(md_files)} markdown files...")

    # Extract all DOIs
    all_dois = []
    for filepath in sorted(md_files):
        dois = extract_dois_from_file(filepath)
        all_dois.extend(dois)

    # Deduplicate
    unique_dois = {}
    for doi, filepath in all_dois:
        if doi not in unique_dois:
            unique_dois[doi] = filepath

    print(f"Found {len(unique_dois)} unique DOIs\n")

    results = []
    failed = []

    for i, (doi, filepath) in enumerate(sorted(unique_dois.items()), 1):
        print(f"[{i}/{len(unique_dois)}] Checking {doi}...")

        result = check_doi(doi)
        result["filepath"] = str(filepath)
        results.append(result)

        if not result["valid"]:
            print(f"  ✗ {result['error']}")
            failed.append((doi, filepath, result["error"]))
        else:
            print(f"  ✓ {result['authors']}: {result['title'][:60]}...")

        # Rate limit
        if i < len(unique_dois):
            time.sleep(1)

    # Write output
    output_file = Path("benchmarks/outputs/doi_check.txt")
    output_file.parent.mkdir(parents=True, exist_ok=True)

    with open(output_file, "w") as f:
        f.write("DOI Validation Results\n")
        f.write("=" * 70 + "\n\n")

        for result in results:
            f.write(f"DOI: {result['doi']}\n")
            f.write(f"File: {result['filepath']}\n")
            if result["valid"]:
                f.write(f"Status: VALID\n")
                f.write(f"Title: {result['title']}\n")
                f.write(f"Authors: {result['authors']}\n")
            else:
                f.write(f"Status: FAILED\n")
                f.write(f"Error: {result['error']}\n")
            f.write("\n")

        f.write("=" * 70 + "\n")
        f.write(f"Total DOIs: {len(unique_dois)}\n")
        f.write(f"Valid: {len(results) - len(failed)}\n")
        f.write(f"Failed: {len(failed)}\n")

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Total DOIs: {len(unique_dois)}")
    print(f"Valid: {len(results) - len(failed)}")
    print(f"Failed: {len(failed)}")
    print(f"\nResults written to {output_file}")

    if failed:
        print("\n❌ Failed DOIs:")
        for doi, filepath, error in failed:
            print(f"  - {doi} (in {filepath}): {error}")
        return 1
    else:
        print("\n✓ All DOIs are valid!")
        return 0


if __name__ == "__main__":
    exit(main())
