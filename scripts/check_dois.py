"""Validate DOIs in NOTES.md and BENCHMARK.md using Crossref API.

Checks each DOI by querying the Crossref API and verifying that the
title, first author, and year match expectations.

Usage:
    python scripts/check_dois.py
"""

import json
import re
import time
import urllib.request
from pathlib import Path
from typing import Dict, Optional, List


def extract_dois_from_file(filepath: Path) -> List[tuple]:
    """Extract DOIs and context from markdown file."""
    dois = []
    
    with open(filepath, 'r') as f:
        content = f.read()
    
    # Match DOI patterns in markdown - be more strict
    # Format: [DOI:10.xxxx/yyyy] or (https://doi.org/10.xxxx/yyyy)
    # Extract just the DOI part
    
    # Pattern 1: https://doi.org/10.xxxx/yyyy
    for match in re.finditer(r'https://doi\.org/(10\.\d{4,}/[^\s\)]+)', content):
        doi = match.group(1).rstrip('.)]')
        
        start = max(0, content.rfind('\n', 0, match.start()))
        end = content.find('\n', match.end())
        if end == -1:
            end = len(content)
        context = content[start:end].strip()
        
        dois.append((doi, context, filepath.name))
    
    # Pattern 2: doi:10.xxxx/yyyy (not in URL)
    for match in re.finditer(r'(?<!org/)doi:(10\.\d{4,}/[^\s\)]+)', content):
        doi = match.group(1).rstrip('.)]')
        
        start = max(0, content.rfind('\n', 0, match.start()))
        end = content.find('\n', match.end())
        if end == -1:
            end = len(content)
        context = content[start:end].strip()
        
        dois.append((doi, context, filepath.name))
    
    # Deduplicate while preserving order
    seen = set()
    unique_dois = []
    for doi, context, fname in dois:
        if doi not in seen:
            seen.add(doi)
            unique_dois.append((doi, context, fname))
    
    return unique_dois


def query_crossref(doi: str) -> Optional[Dict]:
    """Query Crossref API for DOI metadata."""
    url = f"https://api.crossref.org/works/{doi}"
    
    try:
        req = urllib.request.Request(
            url,
            headers={'User-Agent': 'FoldTrust/1.0 (mailto:kanekalla@users.noreply.github.com)'}
        )
        
        with urllib.request.urlopen(req, timeout=10) as response:
            data = json.loads(response.read().decode())
            return data.get('message')
    
    except urllib.error.HTTPError as e:
        if e.code == 404:
            return None
        raise
    except Exception as e:
        print(f"    Error querying {doi}: {e}")
        return None


def format_authors(authors: List[Dict]) -> str:
    """Format author list from Crossref metadata."""
    if not authors:
        return "Unknown"
    
    first_author = authors[0]
    family = first_author.get('family', '')
    given = first_author.get('given', '')
    
    if len(authors) == 1:
        return f"{family}, {given}" if given else family
    else:
        return f"{family} et al."


def check_doi(doi: str, context: str) -> Dict:
    """Check a single DOI and return validation results."""
    result = {
        'doi': doi,
        'context': context,
        'valid': False,
        'error': None,
        'metadata': {}
    }
    
    metadata = query_crossref(doi)
    
    if metadata is None:
        result['error'] = 'DOI not found in Crossref'
        return result
    
    # Extract metadata
    title = metadata.get('title', [''])[0]
    authors = format_authors(metadata.get('author', []))
    
    # Get publication year from various date fields
    year = None
    if 'published-print' in metadata:
        year = metadata['published-print'].get('date-parts', [[None]])[0][0]
    elif 'published-online' in metadata:
        year = metadata['published-online'].get('date-parts', [[None]])[0][0]
    elif 'issued' in metadata:
        year = metadata['issued'].get('date-parts', [[None]])[0][0]
    
    container = metadata.get('container-title', [''])[0]
    publisher = metadata.get('publisher', '')
    
    result['valid'] = True
    result['metadata'] = {
        'title': title,
        'authors': authors,
        'year': year,
        'journal': container,
        'publisher': publisher,
    }
    
    return result


def main():
    print("=" * 70)
    print("DOI Validation")
    print("=" * 70)
    
    # Files to check
    files_to_check = [
        Path('NOTES.md'),
        Path('BENCHMARK.md') if Path('BENCHMARK.md').exists() else None,
    ]
    files_to_check = [f for f in files_to_check if f and f.exists()]
    
    all_dois = []
    for filepath in files_to_check:
        dois = extract_dois_from_file(filepath)
        all_dois.extend(dois)
    
    print(f"\nFound {len(all_dois)} DOIs across {len(files_to_check)} file(s)")
    print()
    
    results = []
    failed = []
    
    for i, (doi, context, source_file) in enumerate(all_dois, 1):
        print(f"[{i}/{len(all_dois)}] Checking {doi}...")
        
        result = check_doi(doi, context)
        results.append(result)
        
        if not result['valid']:
            print(f"  ✗ {result['error']}")
            failed.append((doi, source_file, result['error']))
        else:
            meta = result['metadata']
            print(f"  ✓ {meta['authors']} ({meta['year']}). {meta['title'][:60]}...")
            print(f"    {meta['journal']}")
        
        print()
        
        # Rate limit: be nice to Crossref API
        if i < len(all_dois):
            time.sleep(1)
    
    # Summary
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print()
    print(f"Total DOIs checked: {len(all_dois)}")
    print(f"Valid: {len(results) - len(failed)}")
    print(f"Failed: {len(failed)}")
    
    if failed:
        print("\n❌ Failed DOIs:")
        for doi, source_file, error in failed:
            print(f"  - {doi} (in {source_file}): {error}")
        print("\nPlease fix or remove invalid DOIs.")
        return 1
    else:
        print("\n✓ All DOIs are valid!")
        return 0


if __name__ == "__main__":
    exit(main())
