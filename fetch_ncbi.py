#!/usr/bin/env python3
"""Fetch sequences from NCBI for case rebuilding."""
import urllib.request
import urllib.error
import time
import sys

def fetch_fasta(accession, cache_dir="data/_cache/ncbi"):
    """Fetch FASTA from NCBI using efetch."""
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={accession}&rettype=fasta&retmode=text"
    cache_file = f"{cache_dir}/{accession}.fa"
    
    try:
        print(f"Fetching {accession}...", file=sys.stderr)
        with urllib.request.urlopen(url, timeout=30) as response:
            content = response.read().decode('utf-8')
        
        with open(cache_file, 'w') as f:
            f.write(content)
        
        print(f"Saved to {cache_file}", file=sys.stderr)
        return content
    except urllib.error.URLError as e:
        print(f"Error fetching {accession}: {e}", file=sys.stderr)
        return None

def extract_sequence(fasta_content):
    """Extract just the sequence from FASTA (no header)."""
    lines = fasta_content.strip().split('\n')
    return ''.join(line for line in lines if not line.startswith('>'))

def extract_slice(fasta_content, start, end, strand='+'):
    """Extract slice from FASTA (1-based inclusive coordinates)."""
    seq = extract_sequence(fasta_content)
    # 1-based inclusive to 0-based
    slice_seq = seq[start-1:end]
    
    if strand == '-':
        # Reverse complement
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                     'a': 't', 't': 'a', 'g': 'c', 'c': 'g'}
        slice_seq = ''.join(complement.get(b, b) for b in reversed(slice_seq))
    
    return slice_seq

if __name__ == '__main__':
    cases = [
        ('NG_008728.1', 31999, 32152, '+', 'smn2-iss-n1'),
        ('NM_000492.4', 1, 200, '+', 'cftr-5utr'),
        ('NG_007398.2', 120818, 121000, '+', 'mapt-e10'),
        ('AF009606.1', 44, 118, '+', 'hcv-ires-dii'),
        ('NC_045512.2', 13462, 13542, '+', 'sars2-fse'),
    ]
    
    for accession, start, end, strand, case_name in cases:
        fasta = fetch_fasta(accession)
        if fasta:
            slice_seq = extract_slice(fasta, start, end, strand)
            # Convert to RNA (U)
            rna_seq = slice_seq.replace('T', 'U').replace('t', 'u')
            length = end - start + 1
            
            print(f"\n{case_name}: {accession}:{start}-{end}({strand})")
            print(f"Length: {length} (got {len(slice_seq)})")
            print(f"DNA: {slice_seq[:60]}...")
            print(f"RNA: {rna_seq[:60]}...")
            
            time.sleep(0.5)  # Be nice to NCBI
