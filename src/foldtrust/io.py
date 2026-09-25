"""Input/Output functions for FoldTrust.

Functions for loading sequences, reference structures, and experimental data.
"""

from pathlib import Path
from typing import Optional, Set, Union

import numpy as np
import pandas as pd

from foldtrust._core import FoldData, FoldDataCollection


def parse_dotbracket(structure: str) -> Set[tuple]:
    """
    Parse dot-bracket notation to extract base pairs.

    Supports (), [], {}, <> bracket types for pseudoknots.

    Parameters
    ----------
    structure : str
        Structure in dot-bracket notation

    Returns
    -------
    set of (int, int)
        Set of base pairs (i, j) where i < j (0-indexed)
    """
    pairs = set()
    stacks = {"(": [], "[": [], "{": [], "<": []}
    closers = {")": "(", "]": "[", "}": "{", ">": "<"}

    for i, char in enumerate(structure):
        if char in stacks:
            stacks[char].append(i)
        elif char in closers:
            opener = closers[char]
            if not stacks[opener]:
                raise ValueError(f"Unbalanced bracket at position {i}")
            j = stacks[opener].pop()
            pairs.add((j, i) if j < i else (i, j))

    # Check for unbalanced brackets
    for opener, stack in stacks.items():
        if stack:
            raise ValueError(f"Unbalanced '{opener}' brackets")

    return pairs


def parse_bpseq(filename: Union[str, Path]) -> Set[tuple]:
    """
    Parse bpseq file format.

    Format: each line is "position base pairing_partner"

    Parameters
    ----------
    filename : str or Path
        Path to bpseq file

    Returns
    -------
    set of (int, int)
        Set of base pairs (i, j) where i < j (0-indexed)
    """
    pairs = set()
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 3:
                i = int(parts[0]) - 1  # Convert to 0-indexed
                j = int(parts[2]) - 1
                if j >= 0:  # 0 means unpaired
                    pairs.add((min(i, j), max(i, j)))
    return pairs


def parse_ct(filename: Union[str, Path]) -> Set[tuple]:
    """
    Parse CT (Connectivity Table) file format.

    Format:
    Line 1: N [header]
    Lines 2-N+1: i base i-1 i+1 j i
    where j is the pairing partner (0 if unpaired)

    Parameters
    ----------
    filename : str or Path
        Path to CT file

    Returns
    -------
    set of (int, int)
        Set of base pairs (i, j) where i < j (0-indexed)
    """
    pairs = set()
    with open(filename) as f:
        lines = f.readlines()
        for line in lines[1:]:  # Skip header
            parts = line.strip().split()
            if len(parts) >= 5:
                i = int(parts[0]) - 1  # Convert to 0-indexed
                j = int(parts[4]) - 1
                if j >= 0 and i < j:
                    pairs.add((i, j))
    return pairs


def read_fasta(filename: Union[str, Path], name: Optional[str] = None) -> FoldData:
    """
    Read a FASTA file into a FoldData object.

    Parameters
    ----------
    filename : str or Path
        Path to FASTA file
    name : str, optional
        Name for the FoldData object (defaults to filename stem)

    Returns
    -------
    FoldData
        FoldData object with sequence loaded

    Examples
    --------
    >>> fd = ft.io.read_fasta("data/cases/sars2-fse/sequence.fa")
    >>> fd.sequence[:20]
    'CGGGTTTGCGGTGTAAGTGC'
    """
    filename = Path(filename)

    if name is None:
        name = filename.stem

    with open(filename) as f:
        lines = f.readlines()

    # Skip header lines starting with '>'
    sequence = "".join(line.strip() for line in lines if not line.startswith(">"))

    fd = FoldData(sequence=sequence, name=name)
    fd.uns["source_file"] = str(filename)

    return fd


def read_case_directory(case_dir: Union[str, Path]) -> FoldData:
    """
    Read a disease case directory (sequence.fa + meta.yaml).

    Parameters
    ----------
    case_dir : str or Path
        Path to case directory

    Returns
    -------
    FoldData
        FoldData object with sequence and metadata

    Examples
    --------
    >>> fd = ft.io.read_case_directory("data/cases/sars2-fse")
    >>> fd.uns['disease']
    'COVID-19'
    """
    case_dir = Path(case_dir)
    name = case_dir.name

    # Read sequence
    fd = read_fasta(case_dir / "sequence.fa", name=name)

    # Read metadata if available
    meta_file = case_dir / "meta.yaml"
    if meta_file.exists():
        import yaml

        with open(meta_file) as f:
            metadata = yaml.safe_load(f)
        fd.uns.update(metadata)

    return fd


def read_case_batch(cases_dir: Union[str, Path]) -> FoldDataCollection:
    """
    Read all case directories into a FoldDataCollection.

    Parameters
    ----------
    cases_dir : str or Path
        Path to directory containing case subdirectories

    Returns
    -------
    FoldDataCollection
        Collection of FoldData objects

    Examples
    --------
    >>> collection = ft.io.read_case_batch("data/cases")
    >>> len(collection)
    5
    >>> list(collection.data.keys())
    ['sars2-fse', 'smn2-iss-n1', 'cftr-5utr', 'mapt-e10', 'hcv-ires-dii']
    """
    cases_dir = Path(cases_dir)
    collection = FoldDataCollection()

    for case_dir in sorted(cases_dir.iterdir()):
        if case_dir.is_dir() and (case_dir / "sequence.fa").exists():
            fd = read_case_directory(case_dir)
            collection.add(fd, name=case_dir.name)

    return collection


def read_reference_structure(
    sequence_file: Union[str, Path], structure_file: Union[str, Path], name: Optional[str] = None
) -> FoldData:
    """
    Read a sequence with known reference structure.

    Parameters
    ----------
    sequence_file : str or Path
        Path to FASTA file
    structure_file : str or Path
        Path to structure file (dot-bracket notation)
    name : str, optional
        Name for the FoldData object

    Returns
    -------
    FoldData
        FoldData object with sequence and reference structure

    Examples
    --------
    >>> fd = ft.io.read_reference_structure(
    ...     "benchmarks/reference_data/5S_rRNA_ecoli.fa",
    ...     "benchmarks/reference_data/5S_rRNA_ecoli.ct"
    ... )
    """
    fd = read_fasta(sequence_file, name=name)

    # Read structure (assume dot-bracket in text file)
    with open(structure_file) as f:
        structure = f.read().strip()

    fd.structures["reference"] = structure
    fd.uns["has_reference"] = True

    return fd


def load_shape_data(
    fd: FoldData,
    shape_file: Union[str, Path],
    start: int,
    end: int,
    dataset_name: str = "shape",
) -> FoldData:
    """
    Load SHAPE reactivity data for a genomic region.

    Parameters
    ----------
    fd : FoldData
        FoldData object (sequence must match coordinates)
    shape_file : str or Path
        Path to CSV file with SHAPE data
    start : int
        Start coordinate (1-indexed, inclusive)
    end : int
        End coordinate (1-indexed, inclusive)
    dataset_name : str
        Name for the reactivity column (e.g., 'shape_zhang', 'shape_incarnato')

    Returns
    -------
    FoldData
        Updated FoldData with reactivity in obs

    Examples
    --------
    >>> fd = ft.io.read_fasta("data/cases/sars2-fse/sequence.fa")
    >>> fd = ft.io.load_shape_data(
    ...     fd,
    ...     "/tmp/SARS_CoV-2_shape_comparison/SHAPE data/zhang_invivo_reactivity.csv",
    ...     13468, 13638,
    ...     "shape_zhang"
    ... )
    >>> 'shape_zhang' in fd.obs.columns
    True
    """
    # Read SHAPE data
    shape_df = pd.read_csv(shape_file)

    # Extract reactivities for region
    region_length = end - start + 1
    reactivities = np.full(region_length, np.nan)

    # Detect column names (typically 'Nucleotide' and 'Reactivity' or similar)
    pos_col = [c for c in shape_df.columns if "nucleotide" in c.lower() or "position" in c.lower()][
        0
    ]
    react_col = [c for c in shape_df.columns if "reactivity" in c.lower()][0]

    for i, pos in enumerate(range(start, end + 1)):
        row = shape_df[shape_df[pos_col] == pos]
        if len(row) > 0:
            reactivities[i] = row[react_col].values[0]

    # Add to obs
    fd.obs[dataset_name] = reactivities
    fd.uns[f"{dataset_name}_source"] = str(shape_file)
    fd.uns[f"{dataset_name}_coordinates"] = f"{start}-{end}"

    return fd
