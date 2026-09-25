"""Dataset fetching and parsing for benchmarking."""

import hashlib
import json
import urllib.request
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import gzip


def compute_sha256(filepath: Path) -> str:
    """Compute SHA256 checksum of a file."""
    sha256 = hashlib.sha256()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            sha256.update(chunk)
    return sha256.hexdigest()


def download_file(url: str, dest: Path, expected_sha256: Optional[str] = None) -> None:
    """Download file with optional checksum verification."""
    if dest.exists():
        if expected_sha256 and compute_sha256(dest) == expected_sha256:
            print(f"  ✓ {dest.name} already exists with correct checksum")
            return
        elif not expected_sha256:
            print(f"  ✓ {dest.name} already exists")
            return
    
    print(f"  Downloading {url}")
    dest.parent.mkdir(parents=True, exist_ok=True)
    
    urllib.request.urlretrieve(url, dest)
    
    if expected_sha256:
        actual = compute_sha256(dest)
        if actual != expected_sha256:
            dest.unlink()
            raise ValueError(f"Checksum mismatch for {dest.name}")
        print(f"  ✓ Downloaded and verified {dest.name}")
    else:
        print(f"  ✓ Downloaded {dest.name}")


def fetch_archiveii(data_dir: Path, max_length: Optional[int] = 500) -> Path:
    """
    Fetch ArchiveII dataset (or subset).
    
    ArchiveII is a curated database of RNA structures from comparative analysis.
    Source: RNA STRAND v2.0 / ArchiveII subset
    
    For Mac-16GB friendliness, we use a filtered subset.
    
    Returns:
        Path to the local dataset file
    """
    archive_dir = data_dir / "archiveii"
    archive_dir.mkdir(parents=True, exist_ok=True)
    
    output_file = archive_dir / "archiveii_subset.json"
    
    if output_file.exists():
        print(f"  ✓ ArchiveII subset already exists at {output_file}")
        return output_file
    
    print("  Fetching RNA STRAND database (ArchiveII subset)...")
    
    url = "https://www.rnasoft.ca/strand/RNAStrands.txt"
    raw_file = archive_dir / "RNAStrands.txt"
    
    try:
        download_file(url, raw_file)
    except Exception as e:
        print(f"  Warning: Could not download from rnasoft.ca: {e}")
        print("  Creating minimal test dataset instead...")
        create_minimal_reference_set(output_file)
        return output_file
    
    entries = parse_rna_strand(raw_file, max_length=max_length)
    
    with open(output_file, "w") as f:
        json.dump(entries, f, indent=2)
    
    print(f"  ✓ Parsed {len(entries)} sequences (max length {max_length})")
    
    return output_file


def parse_rna_strand(filepath: Path, max_length: Optional[int] = None) -> List[Dict]:
    """Parse RNA STRAND format file."""
    entries = []
    
    with open(filepath, "r") as f:
        lines = f.readlines()
    
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        
        if line.startswith(">"):
            name = line[1:].strip()
            i += 1
            
            sequence = ""
            structure = ""
            
            while i < len(lines) and not lines[i].startswith(">"):
                content = lines[i].strip()
                if not content:
                    i += 1
                    continue
                
                if all(c in "ACGTU" for c in content.upper()):
                    sequence += content.upper().replace("U", "T")
                elif all(c in ".()[]{}|&" for c in content):
                    structure += content.replace("|", ".").replace("&", ".")
                
                i += 1
            
            if sequence and structure and len(sequence) == len(structure):
                if max_length is None or len(sequence) <= max_length:
                    entries.append({
                        "name": name,
                        "sequence": sequence,
                        "structure": structure,
                        "length": len(sequence),
                        "source": "RNA_STRAND"
                    })
        else:
            i += 1
    
    return entries


def create_minimal_reference_set(output_file: Path) -> None:
    """Create a minimal test reference set for benchmarking."""
    test_sequences = [
        {
            "name": "tRNA_Phe",
            "sequence": "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA",
            "structure": "(((((((..((((.........)))).(((((.......))))).....(((((.......))))))))))))..",
            "length": 76,
            "source": "Manual_test"
        },
        {
            "name": "5S_rRNA_E_coli",
            "sequence": "GGCCUGGCGGCCGUAGCGCGGUGGUCCCACCUGACCCCAUGCCGAACUCAGAAGUGAAACGCCGUAGCGCCGAUGGUAGUGUGGGGUCUCCCCAUGCGAGAGUAGGGAACUGCCAGGCAU",
            "structure": "(((((((((...(((((.......)))))......((((((..(((...)))..))))))......((((((.........))))))......))))))))).................",
            "length": 120,
            "source": "Manual_test"
        },
        {
            "name": "test_hairpin",
            "sequence": "CGCGAAACGCG",
            "structure": "((((...)...))",
            "length": 11,
            "source": "Manual_test"
        }
    ]
    
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, "w") as f:
        json.dump(test_sequences, f, indent=2)
    
    print(f"  ✓ Created minimal test reference set with {len(test_sequences)} sequences")


def load_reference_dataset(filepath: Path) -> List[Dict]:
    """Load reference dataset from JSON."""
    with open(filepath, "r") as f:
        return json.load(f)


def fetch_shape_data_sars2_fse(data_dir: Path) -> Optional[Path]:
    """
    Fetch SHAPE-MaP data for SARS-CoV-2 frameshift element.
    
    Source: Huston et al. 2021, doi:10.1371/journal.ppat.1009265
    or Manfredonia et al. 2020, doi:10.1038/s41586-020-2681-1
    
    Returns path to reactivity file or None if unavailable.
    """
    shape_dir = data_dir / "shape"
    shape_dir.mkdir(parents=True, exist_ok=True)
    
    output_file = shape_dir / "sars2_fse_shape.json"
    
    if output_file.exists():
        print(f"  ✓ SARS2 FSE SHAPE data already exists")
        return output_file
    
    print("  Note: Public SHAPE data for SARS-CoV-2 FSE requires manual extraction from papers")
    print("  Skipping SHAPE analysis for SARS2-FSE (no open single-file source available)")
    
    return None


def create_mock_shape_data(output_file: Path, sequence_length: int) -> None:
    """Create mock SHAPE data for testing (DO NOT USE FOR REAL RESULTS)."""
    import numpy as np
    
    mock_data = {
        "sequence_length": sequence_length,
        "reactivities": list(np.random.rand(sequence_length) * 2),
        "source": "MOCK_DATA_FOR_TESTING_ONLY",
        "note": "This is synthetic data for code testing only"
    }
    
    with open(output_file, "w") as f:
        json.dump(mock_data, f, indent=2)
