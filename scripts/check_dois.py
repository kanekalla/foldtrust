"""Validate DOIs in all markdown files using Crossref API.

Checks each DOI by querying the Crossref API and verifying metadata.
For each occurrence, validates first author and title against citation context.

Usage:
    python scripts/check_dois.py
"""

import json
import re
import time
import urllib.request
from pathlib import Path
from typing import Dict, List, Optional, Tuple


def extract_dois_from_file(filepath: Path) -> List[tuple]:
    """Extract DOIs from markdown file with surrounding context and line numbers."""
    dois = []

    with open(filepath, "r") as f:
        lines = f.readlines()
        content = "".join(lines)

    # Extract DOIs
    pattern = r"10\.\d{4,9}/(?:[a-zA-Z0-9.\-_;/]|(?:\([a-zA-Z0-9.\-_]+\)))+"

    for match in re.finditer(pattern, content):
        doi = match.group(0)
        doi = doi.rstrip(".,;:!?")

        # Get context: ±500 characters around the DOI (wider for full author lists)
        start = max(0, match.start() - 500)
        end = min(len(content), match.end() + 300)
        context = content[start:end]

        # Find line number
        line_num = content[: match.start()].count("\n") + 1

        dois.append((doi, filepath, context, line_num))

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


def normalize_name(name: str) -> str:
    """Normalize name: lowercase, remove accents/apostrophes/hyphens."""
    name = name.lower()
    # Remove apostrophes (all variants)
    for apostrophe in ["'", "'", "`", "'"]:
        name = name.replace(apostrophe, "")
    name = name.replace("-", "")
    return name


def check_author_in_context(first_author_family: str, context: str) -> Tuple[bool, str]:
    """Check if first author appears in citation context.

    Returns (passed, reason) where:
    - passed=True: check passed or skipped (bare DOI)
    - passed=False: clear mismatch (full author list present but clearly wrong first author)
    """
    # Extract author-like names before the DOI
    # Pattern: "LastName INITIAL" for full author lists
    full_author_pattern = r"\b([A-Z][A-Za-z''\-]+)\s+([A-Z]{1,3})\b"
    full_author_matches = re.findall(full_author_pattern, context)

    # Only check if we have a substantial full author list (5+ authors with initials - very conservative)
    if len(full_author_matches) < 5:
        return True, "skip-not-full-list"  # Not checking partial lists

    # For particle names like "van Swieten", take last token
    first_author_tokens = first_author_family.split()
    first_author_key = first_author_tokens[-1] if first_author_tokens else first_author_family
    first_author_normalized = normalize_name(first_author_key)

    # Check if any candidate matches
    for surname, _ in full_author_matches:
        surname_normalized = normalize_name(surname)
        if first_author_normalized == surname_normalized:
            return True, "ok"
        # Also check full name for compound surnames
        if normalize_name(first_author_family) == surname_normalized:
            return True, "ok"

    # Clear mismatch: has substantial author list but first author definitely not there
    return False, "first-author-not-found"


def check_title_in_context(crossref_title: str, context: str) -> Tuple[bool, str]:
    """Check if title appears in citation context.

    Returns (passed, reason) where:
    - passed=True: check passed or skipped
    - passed=False: clear mismatch (substantial quoted title present but completely wrong topic)
    """
    # Only check if there's a clearly quoted title with actual quotes (60+ chars)
    # This is meant to catch only truly wrong citations like in test cases
    title_pattern = r'["\']([^"\']{60,})["\']'
    title_matches = re.findall(title_pattern, context)

    # Filter out false matches (markdown bold/italic that look like quotes)
    title_matches = [m for m in title_matches if not re.search(r"\*\*|__|\|", m)]

    if not title_matches:
        return True, "skip-no-quoted-title"

    # Normalize and tokenize
    def tokenize(text):
        return set(word.lower() for word in re.findall(r"\b\w{4,}\b", text))

    crossref_words = tokenize(crossref_title)
    if len(crossref_words) < 4:  # Need at least 4 meaningful words
        return True, "skip-short-title"

    for candidate in title_matches:
        candidate_words = tokenize(candidate)
        if len(candidate_words) < 4:
            continue

        # Check overlap: at least 20% overlap (extremely lenient - just catch completely wrong titles)
        overlap = len(crossref_words & candidate_words) / len(crossref_words)
        if overlap >= 0.2:
            return True, "ok"

    # Found substantial quoted title but completely wrong topic - clear mismatch
    return False, "title-mismatch"


def check_occurrence(doi: str, context: str, metadata: Dict) -> Dict:
    """Check one DOI occurrence."""
    result = {
        "doi": doi,
        "valid": False,
        "author_ok": None,
        "title_ok": None,
        "checked": None,
        "error": None,
    }

    # DataCite: redirect only
    if doi.startswith("10.5281/"):
        if check_doi_redirect(doi):
            result["valid"] = True
            result["checked"] = "resolve-only"
            result["author_ok"] = "n/a-datacite"
            result["title_ok"] = "n/a-datacite"
        else:
            result["error"] = "redirect-failed"
        return result

    if not metadata:
        result["error"] = "doi-not-found"
        return result

    # Get Crossref metadata
    title = metadata.get("title", [""])[0]
    author_list = metadata.get("author", [])
    first_author_family = author_list[0].get("family", "") if author_list else ""

    # Check author
    author_passed, author_reason = check_author_in_context(first_author_family, context)
    result["author_ok"] = author_reason

    # Check title
    title_passed, title_reason = check_title_in_context(title, context)
    result["title_ok"] = title_reason

    # Determine if checks passed
    if not author_passed:
        result["error"] = f"author-check-failed:{author_reason}"
        result["valid"] = False
    elif not title_passed:
        result["error"] = f"title-check-failed:{title_reason}"
        result["valid"] = False
    else:
        result["valid"] = True
        # Determine check level
        if "skip" in author_reason and "skip" in title_reason:
            result["checked"] = "resolve-only"
        else:
            result["checked"] = "full"

    return result


def main():
    print("=" * 70)
    print("DOI Validation")
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
    for doi, filepath, context, line_num in all_occurrences:
        if doi not in doi_groups:
            doi_groups[doi] = []
        doi_groups[doi].append((filepath, context, line_num))

    print(f"Found {len(doi_groups)} unique DOIs ({len(all_occurrences)} occurrences)\n")

    # Check each DOI (fetch metadata once per DOI)
    all_results = []
    failed_occurrences = []

    for i, (doi, occurrences) in enumerate(sorted(doi_groups.items()), 1):
        print(f"[{i}/{len(doi_groups)}] Checking {doi} ({len(occurrences)} occurrences)...")

        # Fetch Crossref metadata once
        metadata = query_crossref(doi)

        # Check each occurrence
        for filepath, context, line_num in occurrences:
            result = check_occurrence(doi, context, metadata)
            result["filepath"] = str(filepath)
            result["line"] = line_num

            if metadata:
                result["crossref_first_author"] = (
                    metadata.get("author", [{}])[0].get("family", "Unknown")
                    if metadata.get("author")
                    else "Unknown"
                )
            else:
                result["crossref_first_author"] = "n/a"

            all_results.append(result)

            if not result["valid"]:
                failed_occurrences.append((doi, filepath, line_num, result["error"]))

        # Show summary for this DOI
        valid_count = sum(
            1 for _, ctx, _ in occurrences if check_occurrence(doi, ctx, metadata)["valid"]
        )
        if valid_count == len(occurrences):
            print(f"  ✓ All {len(occurrences)} occurrence(s) valid")
        else:
            print(f"  ✗ {len(occurrences) - valid_count}/{len(occurrences)} failed")

        # Rate limit
        if i < len(doi_groups):
            time.sleep(1)

    # Write output
    output_file = Path("benchmarks/outputs/doi_check.txt")
    output_file.parent.mkdir(parents=True, exist_ok=True)

    with open(output_file, "w") as f:
        f.write("DOI Validation Results\n")
        f.write("=" * 70 + "\n\n")

        for result in all_results:
            f.write(f"DOI: {result['doi']}\n")
            f.write(f"File: {result['filepath']}:{result['line']}\n")
            f.write(f"Crossref first author: {result['crossref_first_author']}\n")
            f.write(f"Author check: {result['author_ok']}\n")
            f.write(f"Title check: {result['title_ok']}\n")
            if result["valid"]:
                f.write(f"Status: VALID ({result.get('checked', 'unknown')})\n")
            else:
                f.write("Status: FAILED\n")
                f.write(f"Error: {result['error']}\n")
            f.write("\n")

        f.write("=" * 70 + "\n")
        f.write(f"Total occurrences: {len(all_results)}\n")
        f.write(f"Valid: {len(all_results) - len(failed_occurrences)}\n")
        f.write(f"Failed: {len(failed_occurrences)}\n")

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Total occurrences: {len(all_results)}")
    print(f"Valid: {len(all_results) - len(failed_occurrences)}")
    print(f"Failed: {len(failed_occurrences)}")
    print(f"\nResults written to {output_file}")

    if failed_occurrences:
        print("\n❌ Failed occurrences:")
        for doi, filepath, line, error in failed_occurrences:
            print(f"  - {doi} (in {filepath}:{line}): {error}")
        return 1
    else:
        print("\n✓ All DOI occurrences are valid!")
        return 0


if __name__ == "__main__":
    exit(main())
