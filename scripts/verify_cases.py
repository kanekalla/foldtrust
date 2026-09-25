#!/usr/bin/env python3
"""Verify that case sequences match their NCBI provenance coordinates.

This script:
1. Reads each case's sequence.fa and meta.yaml
2. Fetches the accession from NCBI (or uses cached copy)
3. Extracts the declared slice (start-end, strand)
4. Compares with the case sequence (T/U normalized)
5. Checks that header and meta coordinates agree
6. Reports PASS/FAIL for each case
7. Exits non-zero if any case fails
"""

import argparse
import hashlib
import sys
from pathlib import Path
from typing import Optional, Dict, Any

import yaml

try:
    import urllib.request
    import urllib.error
except ImportError:
    print("ERROR: urllib not available", file=sys.stderr)
    sys.exit(1)


def fetch_fasta(accession: str, cache_dir: Path) -> Optional[str]:
    """Fetch FASTA from NCBI or cache."""
    cache_file = cache_dir / f"{accession}.fa"
    
    if cache_file.exists():
        return cache_file.read_text()
    
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={accession}&rettype=fasta&retmode=text"
    
    try:
        with urllib.request.urlopen(url, timeout=30) as response:
            content = response.read().decode('utf-8')
        
        # Cache it
        cache_dir.mkdir(parents=True, exist_ok=True)
        cache_file.write_text(content)
        
        return content
    except urllib.error.URLError as e:
        print(f"ERROR: Failed to fetch {accession}: {e}", file=sys.stderr)
        return None


def extract_sequence(fasta_content: str) -> str:
    """Extract sequence from FASTA (no header)."""
    lines = fasta_content.strip().split('\n')
    return ''.join(line.strip() for line in lines if not line.startswith('>'))


def normalize_seq(seq: str) -> str:
    """Normalize sequence: uppercase, T<->U interchangeable."""
    return seq.upper().replace('U', 'T')


def reverse_complement(seq: str) -> str:
    """Return reverse complement (DNA alphabet)."""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return ''.join(complement.get(b, b) for b in reversed(seq))


def extract_slice(fasta_content: str, start: int, end: int, strand: str) -> str:
    """Extract slice from FASTA (1-based inclusive coordinates)."""
    seq = extract_sequence(fasta_content)
    # 1-based inclusive -> 0-based
    slice_seq = seq[start - 1:end]
    
    if strand == '-':
        slice_seq = reverse_complement(slice_seq)
    
    return slice_seq


def parse_header_coords(header: str) -> Optional[Dict[str, Any]]:
    """Parse accession:start-end(strand) from FASTA header.
    
    Returns dict with accession, start, end, strand or None if unparseable.
    """
    # Expected format: >casename ACC.VER:start-end(strand) description
    parts = header.strip().lstrip('>').split()
    if len(parts) < 2:
        return None
    
    coord_part = parts[1]
    
    # Parse ACC.VER:start-end(strand)
    if ':' not in coord_part:
        return None
    
    accession, rest = coord_part.split(':', 1)
    
    if '(' not in rest or ')' not in rest:
        return None
    
    range_part, strand_part = rest.split('(', 1)
    strand = strand_part.rstrip(')')
    
    if '-' not in range_part:
        return None
    
    start_str, end_str = range_part.split('-', 1)
    
    try:
        start = int(start_str)
        end = int(end_str)
    except ValueError:
        return None
    
    return {
        'accession': accession,
        'start': start,
        'end': end,
        'strand': strand
    }


def verify_case(case_dir: Path, cache_dir: Path, from_cache: Optional[Path] = None) -> tuple[bool, str]:
    """Verify one case.
    
    Returns (passed, message).
    """
    case_name = case_dir.name
    
    sequence_file = case_dir / "sequence.fa"
    meta_file = case_dir / "meta.yaml"
    
    if not sequence_file.exists():
        return False, f"Missing sequence.fa"
    
    if not meta_file.exists():
        return False, f"Missing meta.yaml"
    
    # Read case sequence
    case_fasta = sequence_file.read_text()
    case_seq = extract_sequence(case_fasta)
    case_seq_norm = normalize_seq(case_seq)
    
    # Parse header coordinates
    header_line = case_fasta.strip().split('\n')[0]
    header_coords = parse_header_coords(header_line)
    
    if not header_coords:
        return False, f"Could not parse header coordinates: {header_line}"
    
    # Read meta
    meta = yaml.safe_load(meta_file.read_text())
    
    if 'provenance' not in meta:
        return False, f"Missing provenance block in meta.yaml"
    
    prov = meta['provenance']
    
    required_keys = ['accession', 'start', 'end', 'strand', 'length']
    missing = [k for k in required_keys if k not in prov]
    if missing:
        return False, f"Missing provenance keys: {missing}"
    
    # Check header and meta agree
    if header_coords['accession'] != prov['accession']:
        return False, f"Header accession {header_coords['accession']} != meta {prov['accession']}"
    
    if header_coords['start'] != prov['start']:
        return False, f"Header start {header_coords['start']} != meta {prov['start']}"
    
    if header_coords['end'] != prov['end']:
        return False, f"Header end {header_coords['end']} != meta {prov['end']}"
    
    if header_coords['strand'] != prov['strand']:
        return False, f"Header strand {header_coords['strand']} != meta {prov['strand']}"
    
    expected_length = prov['end'] - prov['start'] + 1
    if prov['length'] != expected_length:
        return False, f"Meta length {prov['length']} != expected {expected_length}"
    
    if len(case_seq) != expected_length:
        return False, f"Sequence length {len(case_seq)} != expected {expected_length}"
    
    # Fetch reference sequence
    accession = prov['accession']
    start = prov['start']
    end = prov['end']
    strand = prov['strand']
    
    if from_cache:
        # Use cached fixture (already sliced)
        fixture_file = from_cache / f"{accession}_{start}-{end}.fa"
        if not fixture_file.exists():
            return False, f"Fixture not found: {fixture_file}"
        
        ref_fasta = fixture_file.read_text()
        ref_slice = extract_sequence(ref_fasta)  # Already the slice
        ref_slice_norm = normalize_seq(ref_slice)
    else:
        # Fetch from NCBI and extract slice
        ref_fasta = fetch_fasta(accession, cache_dir)
        if ref_fasta is None:
            return False, f"Failed to fetch {accession}"
        
        ref_slice = extract_slice(ref_fasta, start, end, strand)
        ref_slice_norm = normalize_seq(ref_slice)
    
    # Compare
    if case_seq_norm != ref_slice_norm:
        return False, f"Sequence mismatch (normalized)\nCase: {case_seq_norm[:50]}...\n Ref: {ref_slice_norm[:50]}..."
    
    return True, "PASS"


def main():
    parser = argparse.ArgumentParser(description="Verify case provenance")
    parser.add_argument('--from-cache', type=Path, help="Use cached fixtures instead of fetching from NCBI")
    parser.add_argument('--cache-dir', type=Path, default=Path('data/_cache/ncbi'),
                       help="Cache directory for NCBI fetches")
    parser.add_argument('--cases-dir', type=Path, default=Path('data/cases'),
                       help="Directory containing case subdirectories")
    
    args = parser.parse_args()
    
    cases_dir = args.cases_dir
    
    if not cases_dir.exists():
        print(f"ERROR: Cases directory not found: {cases_dir}", file=sys.stderr)
        sys.exit(1)
    
    case_dirs = sorted([d for d in cases_dir.iterdir() if d.is_dir()])
    
    if not case_dirs:
        print(f"ERROR: No case directories found in {cases_dir}", file=sys.stderr)
        sys.exit(1)
    
    print("=" * 80)
    print("FoldTrust Case Provenance Verification")
    print("=" * 80)
    print()
    
    results = []
    
    for case_dir in case_dirs:
        case_name = case_dir.name
        passed, message = verify_case(case_dir, args.cache_dir, args.from_cache)
        
        results.append((case_name, passed, message))
        
        status = "✓ PASS" if passed else "✗ FAIL"
        print(f"{status:10} {case_name:20} {message}")
    
    print()
    print("=" * 80)
    
    passed_count = sum(1 for _, passed, _ in results if passed)
    total_count = len(results)
    
    print(f"Results: {passed_count}/{total_count} passed")
    
    if passed_count < total_count:
        print()
        print("FAILED cases:")
        for case_name, passed, message in results:
            if not passed:
                print(f"  {case_name}: {message}")
        sys.exit(1)
    else:
        print()
        print("All cases verified successfully!")
        sys.exit(0)


if __name__ == '__main__':
    main()
