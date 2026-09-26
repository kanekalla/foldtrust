"""Validate DOIs in all markdown files using Crossref API.

Checks each DOI by querying the Crossref API and verifying:
- First author surname matches for every citation that names an author
- Title similarity for reference-style citations

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
    """Extract DOIs from markdown file with their citation context.

    Context is the specific citation containing the DOI, isolated from neighbouring citations.
    """
    dois = []

    with open(filepath, "r", encoding="utf-8") as f:
        lines = f.readlines()
        content = "".join(lines)

    # Extract DOIs
    pattern = r"10\.\d{4,9}/(?:[a-zA-Z0-9.\-_;/]|(?:\([a-zA-Z0-9.\-_]+\)))+"

    for match in re.finditer(pattern, content):
        doi = match.group(0)
        doi = doi.rstrip(".,;:!?")

        # Skip DOIs that appear in URLs (e.g., https://doi.org/10.xxx)
        # Check if preceded by "://doi.org/" or similar
        prefix_start = max(0, match.start() - 20)
        prefix = content[prefix_start : match.start()]
        if re.search(r"https?://doi\.org/$", prefix):
            continue

        # Find line number
        line_num = content[: match.start()].count("\n") + 1

        # Get context: isolate the specific citation, not neighbouring ones
        # Two cases:
        # 1. List item citations: use the entire list item (starts with "- " or "* ")
        # 2. Inline citations: isolate just the segment for this DOI

        # Check if we're in a list item
        line_start = content.rfind("\n", 0, match.start()) + 1
        line_text = content[line_start : match.start()]
        if re.match(r"^\s*[-*]\s+", line_text):
            # List item citation - use the entire item as context
            # Expand backwards if it's a multi-line list item
            item_start = line_start
            while item_start > 0:
                prev_line_end = item_start - 1
                prev_line_start = content.rfind("\n", 0, prev_line_end) + 1
                prev_line = content[prev_line_start:prev_line_end]
                # Continue if the previous line is a continuation (starts with spaces/tabs, not a new list item)
                if prev_line and not re.match(r"^\s*[-*]\s+", prev_line) and prev_line[0] in " \t":
                    item_start = prev_line_start
                else:
                    break

            # Expand forward to end of item
            line_end = content.find("\n", match.end())
            if line_end == -1:
                line_end = len(content)
            item_end = line_end

            # Use the entire list item as context
            context = content[item_start:item_end].strip()
            dois.append((doi, filepath, context, line_num))
            continue

        # Not a list item - extract inline citation context
        context_start = match.start()
        search_limit = max(0, match.start() - 200)

        # Skip back over "doi:" prefix, markdown link syntax, and any preceding semicolon+space
        doi_prefix_start = match.start()

        # Check for markdown link: [10.xxx](url) or DOI [10.xxx](url)
        if match.start() >= 1 and content[match.start() - 1] == "[":
            doi_prefix_start = match.start() - 1
            # Also check for "DOI " before the bracket
            if (
                doi_prefix_start >= 4
                and content[doi_prefix_start - 4 : doi_prefix_start].upper() == "DOI "
            ):
                doi_prefix_start -= 4
        # Check for inline "doi:" prefix
        elif match.start() >= 4 and content[match.start() - 4 : match.start()].lower() == "doi:":
            doi_prefix_start = match.start() - 4
            # Also skip preceding "; " if present (e.g., "; doi:X")
            if doi_prefix_start >= 2 and content[doi_prefix_start - 2 : doi_prefix_start] == "; ":
                doi_prefix_start -= 2
            elif doi_prefix_start >= 1 and content[doi_prefix_start - 1] == ";":
                doi_prefix_start -= 1
        elif match.start() >= 5 and content[match.start() - 5 : match.start()].lower() == " doi:":
            doi_prefix_start = match.start() - 5
            # Also skip preceding ";" if present
            if doi_prefix_start >= 1 and content[doi_prefix_start - 1] == ";":
                doi_prefix_start -= 1

        # Search backward for citation boundary
        # Stop at: semicolon separating citations, opening paren, or list/paragraph boundary
        for i in range(doi_prefix_start - 1, search_limit, -1):
            if content[i] == ";":
                # Semicolon - this separates two citations in the same group
                # The next citation should start with author name (capital letter)
                next_content = content[i + 1 : min(len(content), i + 30)].lstrip()
                if next_content and next_content[0].isupper():
                    context_start = i + 1
                    break
            elif content[i] == "(":
                # Opening paren starts a citation group
                context_start = i
                break
            elif i > 0 and content[i - 1 : i + 1] == ")(":
                # Closing+opening paren indicates previous citation ended
                context_start = i
                break
            elif content[i] == ".":
                # Period - check if it's a sentence boundary (not "et al." or initial)
                if i >= 3 and content[i - 3 : i] == " al":
                    continue
                elif (
                    i >= 1
                    and content[i - 1].isupper()
                    and i + 1 < len(content)
                    and content[i + 1] == " "
                ):
                    continue
                else:
                    # Likely sentence boundary
                    context_start = i + 1
                    break
            elif content[i] == "\n":
                # Newline - check context
                line_start = content.rfind("\n", 0, i) + 1
                prev_line = content[line_start:i]
                if re.match(r"^\s*[-*]\s+", prev_line):
                    context_start = line_start
                    break
                elif not prev_line.strip():
                    context_start = i + 1
                    break
        else:
            context_start = search_limit

        # Look forward for end: period, semicolon, closing paren, or newline
        context_end = match.end()
        for i in range(match.end(), min(len(content), match.end() + 50)):
            if content[i] in ").\n;":
                context_end = i + 1
                break
        else:
            context_end = min(len(content), match.end() + 50)

        context = content[context_start:context_end].strip()
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
    """Normalize name: NFKD, strip accents, lowercase, remove apostrophes/hyphens/spaces."""
    # NFKD normalization and strip accents
    nfkd = unicodedata.normalize("NFKD", name)
    name = "".join(c for c in nfkd if not unicodedata.combining(c))

    name = name.lower()

    # Remove all apostrophe variants, hyphens, backticks, spaces
    for char in ["'", "'", "'", "`", "-", " "]:
        name = name.replace(char, "")

    return name


def extract_first_cited_author(context: str) -> Optional[str]:
    """Extract the first author surname from the citation context."""
    # Try list-style first: "Surname INIT, ..." or "Surname INIT."
    # Allow Unicode apostrophes, hyphens, and spaces in surnames
    list_match = re.search(r"[-*]?\s*([A-Z][\w''\-\u2019 ]+?)\s+[A-Z]{1,3}[,.]", context)
    if list_match:
        surname = list_match.group(1).strip()
        # Filter out common non-name words
        if surname.lower() not in ["nature", "science", "cell", "proc", "natl", "acad"]:
            return surname

    # Try prose style: "Surname et al." or "Surname and ..." BEFORE author-year patterns
    # This prevents matching journal names like "Nature 1998"
    prose_match = re.search(r"\b([A-Z][\w''\-\u2019]+)\s+(?:et\s+al\.?|and\s+[A-Z])", context)
    if prose_match:
        return prose_match.group(1)

    # Try prose style: "Surname, Surname & Surname YEAR" or "Surname, Surname, Surname YEAR"
    # This catches author lists without initials, like "Langdon, Petke & Lorenz 2018"
    multi_author_match = re.search(
        r"\b([A-Z][\w''\-\u2019]+)(?:,\s+[A-Z][\w''\-\u2019]+)+(?:\s+[&,]\s+[A-Z][\w''\-\u2019]+)?\s+(?:19|20)\d{2}\b",
        context,
    )
    if multi_author_match:
        return multi_author_match.group(1)

    # Try standalone author-year: "Surname YEAR"
    year_match = re.search(r"\b([A-Z][\w''\-\u2019]+)\s+(?:19|20)\d{2}\b", context)
    if year_match:
        surname = year_match.group(1)
        if surname.lower() not in ["nature", "science", "cell", "proc", "natl", "acad"]:
            return surname

    return None


def check_author_match(first_author_family: str, context: str) -> Tuple[bool, str]:
    """Check if first cited author matches Crossref first author.

    Returns (passed, reason)
    """
    cited_author = extract_first_cited_author(context)

    if not cited_author:
        return True, "skip-no-author"

    # Normalize both names
    cited_normalized = normalize_name(cited_author)
    crossref_normalized = normalize_name(first_author_family)

    # Direct match
    if cited_normalized == crossref_normalized:
        return True, "ok"

    # For particle names (van Swieten, de Graaff, etc.), check last token
    crossref_tokens = first_author_family.split()
    if len(crossref_tokens) > 1:
        last_token_normalized = normalize_name(crossref_tokens[-1])
        if cited_normalized == last_token_normalized:
            return True, "ok-particle"

    return False, "first-author-mismatch"


def extract_reference_title(context: str) -> Optional[str]:
    """Extract reference-style title from citation context.

    Only extracts titles that are explicitly delimited with quotes or italics.
    Does NOT extract unquoted titles or single-word journal names in italics.
    """
    # Look for quoted title segments
    # Do NOT treat single apostrophe (') as a quote delimiter (to avoid 5' false matches)
    quoted_patterns = [
        r'["""]([^"""]{10,})["""]',  # Double quotes (straight or curly)
    ]

    for pattern in quoted_patterns:
        match = re.search(pattern, context)
        if match:
            return match.group(1).strip()

    # Look for italicized/emphasized title segments (using single * or _)
    # But only if they contain multiple words (to exclude journal names like "*Nature*")
    # Use negative lookbehind/lookahead to avoid matching **bold** (double asterisks)
    emphasized_patterns = [
        r"(?<!\*)\*([^*]{10,})\*(?!\*)",  # *emphasis* but not **bold**
        r"(?<!_)_([^_]{10,})_(?!_)",  # _emphasis_ but not __bold__
    ]

    for pattern in emphasized_patterns:
        match = re.search(pattern, context)
        if match:
            potential_title = match.group(1).strip()
            # Only consider it a title if it has at least 2 words AND doesn't end with a colon (labels)
            if potential_title.endswith(":"):
                continue
            words = potential_title.split()
            if len(words) < 2:
                continue

            # Check if this looks like a journal name (appears before year/volume like "(2011)" or "23(4)")
            # Journal pattern: italicized text followed by optional space and year/volume
            after_match = context[match.end() : match.end() + 20]
            if re.match(r"\s*\(?\d{4}\)?|\s+\d+\(", after_match):
                # Likely a journal name, not a title
                continue

            return potential_title

    return None


def check_title_match(crossref_title: str, context: str) -> Tuple[bool, str]:
    """Check if title appears in citation context.

    Returns (passed, reason)
    """
    cited_title = extract_reference_title(context)

    if not cited_title:
        return True, "skip-no-title"

    # Normalize and tokenize (words of length >= 4)
    def tokenize(text):
        return set(word.lower() for word in re.findall(r"\b\w{4,}\b", text))

    crossref_words = tokenize(crossref_title)
    cited_words = tokenize(cited_title)

    if len(crossref_words) == 0:
        return True, "skip-short-title"

    # Check overlap: at least 50% of Crossref words appear
    overlap = len(crossref_words & cited_words)
    ratio = overlap / len(crossref_words)

    if ratio >= 0.5:
        return True, "ok"

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
    author_passed, author_reason = check_author_match(first_author_family, context)
    result["author_ok"] = author_reason

    # Check title
    title_passed, title_reason = check_title_match(title, context)
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
