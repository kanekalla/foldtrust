"""Utility functions for FoldTrust."""

from pathlib import Path
from typing import List


def read_fasta(fasta_path: Path) -> tuple[str, str]:
    """
    Read sequence from FASTA file.
    
    Returns:
        (header, sequence) tuple
    """
    with open(fasta_path, 'r') as f:
        lines = f.readlines()
    
    header = ""
    sequence_lines = []
    
    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            header = line[1:].strip()
        else:
            sequence_lines.append(line.upper())
    
    sequence = "".join(sequence_lines)
    
    sequence = sequence.replace("U", "T")
    
    return header, sequence


def compute_gc_content(sequence: str) -> float:
    """Calculate GC content as percentage."""
    gc_count = sequence.count("G") + sequence.count("C")
    return (gc_count / len(sequence)) * 100 if sequence else 0.0


def find_case_directories(cases_dir: Path) -> List[Path]:
    """Find all case directories containing sequence.fa."""
    case_dirs = []
    for item in sorted(cases_dir.iterdir()):
        if item.is_dir() and (item / "sequence.fa").exists():
            case_dirs.append(item)
    return case_dirs


def format_structure_ascii(sequence: str, structure: str, width: int = 80) -> str:
    """Format structure in ASCII for display."""
    lines = []
    for i in range(0, len(sequence), width):
        seq_chunk = sequence[i:i+width]
        struct_chunk = structure[i:i+width]
        pos_start = i + 1
        lines.append(f"{pos_start:6d} {seq_chunk}")
        lines.append(f"       {struct_chunk}\n")
    return "\n".join(lines)
