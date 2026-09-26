"""Validate DOIs in all markdown files using Crossref API.

Checks each DOI by querying the Crossref API and verifying:
- DOI resolution for every occurrence (Crossref or doi.org redirect for DataCite)
- First author surname matches for every citation whose segment names an author
- Title similarity only when double-quoted with 4+ words

Reference-style titles (unquoted) are not checked automatically. All titles and
first authors were checked manually (see docs/benchmark/impact.md:157).

Usage:
    python scripts/check_dois.py
"""

import json
import re
import time
import unicodedata
import urllib.request
from pathlib import Path
from typing import Dict, List, Optional, Tuple


def extract_dois_from_file(filepath: Path) -> List[tuple]:
    """Extract DOIs from markdown file with per-citation context segments.

    For each DOI, extracts the citation segment: text on the same line from the
    end of the previous DOI (or line start) to this DOI. For markdown links
    [text](url), uses the link text as the segment.
    """
    dois = []

    with open(filepath, "r", encoding="utf-8") as f:
        lines = f.readlines()
        content = "".join(lines)

    # Extract DOIs
    pattern = r"10\.\d{4,9}/(?:[a-zA-Z0-9.\-_;/]|(?:\([a-zA-Z0-9.\-_]+\)))+"

    # Track DOIs per line for segmentation
    line_dois = {}
    for match in re.finditer(pattern, content):
        doi = match.group(0)
        doi = doi.rstrip(".,;:!?")
        line_num = content[: match.start()].count("\n") + 1

        if line_num not in line_dois:
            line_dois[line_num] = []
        line_dois[line_num].append((doi, match.start(), match.end()))

    # Process each line's DOIs to build segments
    for line_num, line_doi_list in sorted(line_dois.items()):
        line_start = content.rfind("\n", 0, line_doi_list[0][1]) + 1
        line_end = content.find("\n", line_doi_list[-1][2])
        if line_end == -1:
            line_end = len(content)
        line_text = content[line_start:line_end]

        prev_end = 0  # End position relative to line start
        for doi, abs_start, abs_end in line_doi_list:
            # Calculate positions relative to line start
            rel_start = abs_start - line_start
            rel_end = abs_end - line_start

            # R1(i): Check if DOI is in a markdown link [text](url)
            # Look for pattern [text](https://doi.org/DOI) or [text](DOI)
            link_match = None
            for m in re.finditer(r"\[([^\]]+)\]\(([^)]+)\)", line_text):
                link_text, link_url = m.groups()
                # Check if this DOI appears in the URL
                if doi in link_url:
                    link_match = (link_text, m.start(), m.end())
                    break

            if link_match:
                segment, _, _ = link_match
                # Check if link text is itself a DOI
                if re.match(r"(?:DOI:|doi:)?\s*10\.\d", segment):
                    # Duplicate, skip
                    continue
            else:
                # Regular segment: from prev_end to start of "doi:" or DOI number
                segment_start = prev_end
                segment_end = rel_start

                # Check if there's "doi:" or "DOI:" before the number
                prefix_search = line_text[max(0, segment_end - 10) : segment_end]
                if re.search(r"doi:\s*$", prefix_search, re.IGNORECASE):
                    segment_end = max(
                        0,
                        segment_end
                        - len(re.search(r"doi:\s*$", prefix_search, re.IGNORECASE).group()),
                    )

                segment = line_text[segment_start:segment_end]

            # R1(ii): Strip leading ), ], (https...), ;, ., whitespace, bullet
            segment = segment.strip()
            segment = re.sub(r"^[)\];.\s]+", "", segment)
            segment = re.sub(r"^\(https?://[^)]+\)\s*", "", segment)
            segment = re.sub(r"^[-*]\s+", "", segment)

            # R1(iii): If segment contains (, keep only text after last (
            if "(" in segment:
                segment = segment[segment.rfind("(") + 1 :]

            # R1(iv): Strip leading label like "SMN2 ISS-N1: " or "**ViennaRNA**: "
            label_match = re.match(r"^[^:()]{1,40}:\s+(.+)", segment)
            if label_match and re.match(r"^[A-Z]", label_match.group(1)):
                segment = label_match.group(1)

            # R1(v): Remove markdown *, _, **, [, ]
            segment = re.sub(r"\*+|_+|\[|\]", "", segment)
            segment = segment.strip()

            dois.append((doi, filepath, segment, line_num))

            # Update prev_end for next segment
            prev_end = rel_end

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
    """Normalize name: NFKD, strip accents, lowercase, remove apostrophes/hyphens/spaces.

    R4: Remove all apostrophe variants, hyphens, backticks, spaces.
    """
    # NFKD normalization and strip accents
    nfkd = unicodedata.normalize("NFKD", name)
    name = "".join(c for c in nfkd if not unicodedata.combining(c))

    name = name.lower()

    # R4: Remove all apostrophe variants (U+0027, U+2019, U+2018, U+02BC), hyphens, backticks, spaces
    for char in ["'", "\u2019", "\u2018", "\u02bc", "`", "-", " "]:
        name = name.replace(char, "")

    return name


def extract_first_cited_author(segment: str) -> Optional[str]:
    """Extract the first author surname from the citation segment.

    R2: Anchored at segment start, try patterns in order.
    Name token N = [A-Z][\\w'''ʼ\\-]+, optionally preceded by particles.
    """
    # R2: Define name token with optional particles
    # Particles: van, de, der, den, von, zu, "Höner zu"
    name_token = r"(?:(?:van|de|der|den|von|zu|Höner zu)\s+)?[A-Z][\w'''ʼ\-]+"

    # R2.1: "Name INIT" (1-3 uppercase letters)
    match = re.match(rf"^({name_token})\s+[A-Z]{{1,3}}\b", segment)
    if match:
        name = match.group(1)
        # Never return all-caps 1-3 letter tokens
        if not (name.isupper() and len(name) <= 3):
            return name

    # R2.2: "Name et al." or "Name and" or "Name &"
    match = re.match(rf"^({name_token})\s+(?:et\s+al\.?|and\b|&)", segment)
    if match:
        name = match.group(1)
        if not (name.isupper() and len(name) <= 3):
            return name

    # R2.3: "Name, Name, ... & Name" (multi-author, comma-separated)
    match = re.match(
        rf"^({name_token})(?:,\s+{name_token})*\s*(?:,|&|and)\s*{name_token}",
        segment,
    )
    if match:
        name = match.group(1)
        if not (name.isupper() and len(name) <= 3):
            return name

    # R2.4: "Name YEAR"
    match = re.match(rf"^({name_token})\s+(?:19|20)\d{{2}}\b", segment)
    if match:
        name = match.group(1)
        if not (name.isupper() and len(name) <= 3):
            return name

    return None


def check_author_match(first_author_family: str, segment: str) -> Tuple[bool, str]:
    """Check if first cited author matches Crossref first author.

    Returns (passed, reason)
    """
    cited_author = extract_first_cited_author(segment)

    if not cited_author:
        return True, "skip-no-author"

    # Normalize both names
    cited_normalized = normalize_name(cited_author)
    crossref_normalized = normalize_name(first_author_family)

    # Direct match
    if cited_normalized == crossref_normalized:
        return True, "ok"

    # R5: For particle names, check last token
    crossref_tokens = first_author_family.split()
    if len(crossref_tokens) > 1:
        last_token_normalized = normalize_name(crossref_tokens[-1])
        if cited_normalized == last_token_normalized:
            return True, "ok-particle"

    return False, "first-author-mismatch"


def extract_title(segment: str) -> Optional[str]:
    """Extract double-quoted title with 4+ words.

    R6: Only check titles in double quotes ("..." or "...") with 4+ words.
    Never treat italics or single quotes as title delimiters.
    """
    # R6: Look for double quotes only (straight or curly)
    patterns = [
        r'["""]([^"""]+)["""]',  # Double quotes
    ]

    for pattern in patterns:
        match = re.search(pattern, segment)
        if match:
            potential_title = match.group(1).strip()
            # Only consider titles with 4+ words
            words = potential_title.split()
            if len(words) >= 4:
                return potential_title

    return None


def check_title_match(crossref_title: str, segment: str) -> Tuple[bool, str]:
    """Check if title appears in citation segment.

    Returns (passed, reason)
    """
    cited_title = extract_title(segment)

    if not cited_title:
        return True, "skip-no-title"

    # Normalize and tokenize (words of length >= 4)
    def tokenize(text):
        return set(word.lower() for word in re.findall(r"\b\w{4,}\b", text))

    crossref_words = tokenize(crossref_title)
    cited_words = tokenize(cited_title)

    if len(crossref_words) == 0:
        return True, "skip-short-title"

    # R6: Check overlap: at least 50% of Crossref words appear
    overlap = len(crossref_words & cited_words)
    ratio = overlap / len(crossref_words)

    if ratio >= 0.5:
        return True, "ok"

    return False, "title-mismatch"


def check_occurrence(doi: str, segment: str, metadata: Dict) -> Dict:
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
    author_passed, author_reason = check_author_match(first_author_family, segment)
    result["author_ok"] = author_reason

    # Check title
    title_passed, title_reason = check_title_match(title, segment)
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
        if author_reason.startswith("skip") and title_reason.startswith("skip"):
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
    for doi, filepath, segment, line_num in all_occurrences:
        if doi not in doi_groups:
            doi_groups[doi] = []
        doi_groups[doi].append((filepath, segment, line_num))

    print(f"Found {len(doi_groups)} unique DOIs ({len(all_occurrences)} occurrences)\n")

    # Check each DOI (fetch metadata once per DOI)
    all_results = []
    failed_occurrences = []

    for i, (doi, occurrences) in enumerate(sorted(doi_groups.items()), 1):
        print(f"[{i}/{len(doi_groups)}] Checking {doi} ({len(occurrences)} occurrences)...")

        # Fetch Crossref metadata once
        metadata = query_crossref(doi)

        # Check each occurrence
        for filepath, segment, line_num in occurrences:
            result = check_occurrence(doi, segment, metadata)
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
            1 for _, seg, _ in occurrences if check_occurrence(doi, seg, metadata)["valid"]
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
