"""Fetch Rfam seed alignments and extract reference structures.

This script builds a curated reference dataset from Rfam seed alignments,
projecting consensus SS_cons structures onto individual seed sequences.

Rfam families included:
- RF00005: tRNA
- RF00001: 5S rRNA
- RF00002: 5.8S rRNA  
- RF00167: Purine riboswitch
- RF00059: TPP riboswitch (THI element)
- RF00162: SAM riboswitch (S box)
- RF00050: FMN riboswitch
- RF00017: SRP_bact (Bacterial SRP RNA)
- RF01854: SRP_euk_arch (Eukaryotic/Archaeal SRP RNA)
- RF00010: RNase_P_bact_A (Bacterial RNase P class A)
- RF00373: RNase_MRP (RNase MRP)
- RF00012: U3 snoRNA
- RF00015: U4 snRNA
- RF01846: U1 snRNA
- RF00004: U2 snRNA
- RF00026: U6 snRNA

Usage:
    python scripts/build_rfam_dataset.py --output data/rfam_references.json
"""

import argparse
import gzip
import hashlib
import json
import re
import urllib.request
from pathlib import Path
from typing import Dict, List, Optional, Tuple


# Rfam families to include (family_id, max_sequences, max_length)
RFAM_FAMILIES = [
    ("RF00005", "tRNA", 20, 200),
    ("RF00001", "5S_rRNA", 15, 200),
    ("RF00002", "5.8S_rRNA", 10, 200),
    ("RF00167", "Purine_riboswitch", 10, 150),
    ("RF00059", "TPP_riboswitch", 10, 150),
    ("RF00162", "SAM_riboswitch", 10, 200),
    ("RF00050", "FMN_riboswitch", 10, 200),
    ("RF00017", "SRP_bact", 10, 300),
    ("RF01854", "SRP_euk_arch", 5, 350),
    ("RF00010", "RNase_P_bact_A", 5, 400),
    ("RF00012", "U3_snoRNA", 5, 300),
]


def compute_sha256(data: bytes) -> str:
    """Compute SHA256 checksum."""
    return hashlib.sha256(data).hexdigest()


def download_with_checksum(url: str, output_path: Path) -> str:
    """Download file and return its SHA256 checksum."""
    print(f"  Downloading {url}")
    response = urllib.request.urlopen(url)
    data = response.read()
    
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'wb') as f:
        f.write(data)
    
    checksum = compute_sha256(data)
    print(f"  SHA256: {checksum}")
    return checksum


def parse_stockholm(content: str) -> Optional[Tuple[str, str, List[Tuple[str, str, str]]]]:
    """
    Parse Stockholm format seed alignment.
    
    Returns:
        (family_id, family_name, [(seq_id, sequence, structure), ...])
    """
    lines = content.strip().split('\n')
    
    family_id = None
    family_name = None
    ss_cons = None
    sequences = []
    seq_dict = {}
    
    for line in lines:
        line = line.rstrip()
        
        if line.startswith('#=GF ID'):
            family_name = line.split(None, 2)[2]
        elif line.startswith('#=GF AC'):
            family_id = line.split(None, 2)[2].split('.')[0]  # Remove version
        elif line.startswith('#=GC SS_cons'):
            # Consensus structure
            parts = line.split(None, 2)
            if len(parts) >= 3:
                ss_cons = parts[2]
        elif line.startswith('#') or line.startswith('//') or not line:
            continue
        else:
            # Sequence line
            parts = line.split(None, 1)
            if len(parts) == 2:
                seq_id, seq_with_gaps = parts
                if seq_id not in seq_dict:
                    seq_dict[seq_id] = ""
                seq_dict[seq_id] += seq_with_gaps
    
    if not ss_cons or not seq_dict:
        return None
    
    # Project SS_cons onto each sequence (remove gaps)
    for seq_id, seq_with_gaps in seq_dict.items():
        ungapped_seq = seq_with_gaps.replace('.', '').replace('-', '').upper().replace('U', 'T')
        
        # Project structure (remove positions where sequence has gaps)
        projected_structure = ""
        seq_pos = 0
        for i, char in enumerate(seq_with_gaps):
            if char not in '.—':
                if i < len(ss_cons):
                    ss_char = ss_cons[i]
                    # Convert structure notation
                    if ss_char in '<':
                        ss_char = '('
                    elif ss_char in '>':
                        ss_char = ')'
                    elif ss_char in ':,_-~':
                        ss_char = '.'
                    projected_structure += ss_char
                else:
                    projected_structure += '.'
                seq_pos += 1
        
        # Remove pseudoknot notation (keep only ().)
        clean_structure = ""
        paren_stack = []
        for char in projected_structure:
            if char in '([{<':
                clean_structure += '('
                paren_stack.append('(')
            elif char in ')]}>{':
                if paren_stack:
                    paren_stack.pop()
                    clean_structure += ')'
                else:
                    clean_structure += '.'
            else:
                clean_structure += '.'
        
        if len(ungapped_seq) == len(clean_structure) and ungapped_seq:
            sequences.append((seq_id, ungapped_seq, clean_structure))
    
    return (family_id, family_name, sequences)


def fetch_rfam_family(family_id: str, cache_dir: Path) -> Optional[Tuple[str, List[Dict]]]:
    """
    Fetch Rfam seed alignment for a family.
    
    Returns:
        (sha256, [sequence_dicts])
    """
    # Use Rfam FTP: ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/
    # Format: RF00005.seed.gz
    url = f"https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/{family_id}.seed.gz"
    
    cache_file = cache_dir / f"{family_id}.seed.gz"
    
    if not cache_file.exists():
        try:
            checksum = download_with_checksum(url, cache_file)
        except Exception as e:
            print(f"  Error downloading {family_id}: {e}")
            return None
    else:
        with open(cache_file, 'rb') as f:
            data = f.read()
        checksum = compute_sha256(data)
        print(f"  Using cached {family_id} (SHA256: {checksum})")
    
    # Parse the Stockholm file
    try:
        with gzip.open(cache_file, 'rt') as f:
            content = f.read()
    except Exception as e:
        print(f"  Error reading {cache_file}: {e}")
        return None
    
    result = parse_stockholm(content)
    if not result:
        print(f"  Failed to parse {family_id}")
        return None
    
    family_id_parsed, family_name, sequences = result
    
    entries = []
    for seq_id, sequence, structure in sequences:
        entries.append({
            "name": f"{family_name}_{seq_id}",
            "sequence": sequence,
            "structure": structure,
            "length": len(sequence),
            "source": f"Rfam_{family_id}",
            "rfam_id": family_id,
            "rfam_name": family_name,
        })
    
    return (checksum, entries)


def build_rfam_dataset(output_file: Path, cache_dir: Path, max_total: int = 300):
    """Build reference dataset from Rfam seed alignments."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    
    all_entries = []
    checksums = {}
    
    print(f"\nFetching Rfam families...")
    for family_id, family_name, max_seqs, max_len in RFAM_FAMILIES:
        print(f"\n{family_id} ({family_name}):")
        result = fetch_rfam_family(family_id, cache_dir)
        
        if result:
            checksum, entries = result
            checksums[family_id] = checksum
            
            # Filter by length and sample
            filtered = [e for e in entries if e['length'] <= max_len]
            sampled = filtered[:max_seqs]
            
            print(f"  Collected {len(sampled)}/{len(entries)} sequences (max_len={max_len})")
            all_entries.extend(sampled)
        
        if len(all_entries) >= max_total:
            break
    
    # Save dataset
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, 'w') as f:
        json.dump(all_entries, f, indent=2)
    
    print(f"\n✓ Saved {len(all_entries)} sequences to {output_file}")
    
    # Save checksums
    checksum_file = output_file.parent / "rfam_checksums.json"
    with open(checksum_file, 'w') as f:
        json.dump({
            "source": "Rfam seed alignments from https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/",
            "date": "2026-09-25",
            "families": checksums,
        }, f, indent=2)
    
    print(f"✓ Saved checksums to {checksum_file}")
    
    return all_entries


def main():
    parser = argparse.ArgumentParser(description="Build Rfam reference dataset")
    parser.add_argument("--output", type=Path, default=Path("data/rfam_references.json"))
    parser.add_argument("--cache", type=Path, default=Path("data/rfam_cache"))
    parser.add_argument("--max-total", type=int, default=300)
    
    args = parser.parse_args()
    
    build_rfam_dataset(args.output, args.cache, args.max_total)
    print("\n✓ Done!")


if __name__ == "__main__":
    main()
