"""Layer 5: Robustness analysis across parameter sets, temperatures, and window jitter.

This module implements comprehensive robustness testing of FoldTrust reliability tiers
under three types of perturbations:
1. Energy parameter sets (Turner2004, Andronescu2007, Langdon2018)
2. Temperature variations (24°C, 30°C, 37°C baseline, 42°C, 45°C)
3. Window boundary jitter (extend/shrink by 10 and 25 nt on each side)

All conditions are compared to the baseline: Turner2004, 37°C, original window.
"""

import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import urllib.request
import urllib.error

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from foldtrust.vienna import (
    fold_with_params,
    compute_pair_probs_with_params,
    parse_stems,
    ViennaRNAError,
)
from foldtrust.utils import read_fasta


# Metadata for fetching flanking sequences from NCBI
CASE_GENOMIC_COORDS = {
    "sars2-fse": {
        "accession": "NC_045512.2",
        "start": 13462,  # Slippery site start
        "end": 13542,    # End of stem-loop 2
        "strand": "+",
    },
    "smn2-iss-n1": {
        "accession": "NM_017411.4",
        "start": 840,
        "end": 1040,
        "strand": "+",
    },
    "cftr-5utr": {
        "accession": "NM_000492.4",
        "start": 133,
        "end": 400,
        "strand": "+",
    },
    "mapt-e10": {
        "accession": "NM_001123066.4",
        "start": 980,
        "end": 1260,
        "strand": "+",
    },
    "hcv-ires-dii": {
        "accession": "AF009606.1",
        "start": 40,
        "end": 310,
        "strand": "+",
    },
}


def fetch_genomic_region(
    accession: str,
    start: int,
    end: int,
    strand: str = "+",
) -> str:
    """
    Fetch genomic sequence from NCBI E-utilities.
    
    Args:
        accession: GenBank/RefSeq accession
        start: 1-based start position
        end: 1-based end position (inclusive)
        strand: "+" or "-"
    
    Returns:
        Sequence string
        
    Raises:
        RuntimeError: If fetch fails
    """
    url = (
        f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
        f"db=nuccore&id={accession}&rettype=fasta&retmode=text"
        f"&seq_start={start}&seq_stop={end}"
    )
    
    if strand == "-":
        url += "&strand=2"
    
    try:
        with urllib.request.urlopen(url, timeout=30) as response:
            content = response.read().decode('utf-8')
            lines = content.strip().split('\n')
            # Skip header line
            sequence = ''.join(lines[1:]).upper().replace('U', 'T')
            return sequence.replace('T', 'U')  # Convert to RNA
    except urllib.error.URLError as e:
        raise RuntimeError(f"Failed to fetch {accession}:{start}-{end}: {e}")


def verify_case_sequence(case_name: str, case_sequence: str) -> Tuple[bool, str]:
    """
    Verify that case sequence matches genomic coordinates.
    
    Returns:
        (matches, message) tuple
    """
    if case_name not in CASE_GENOMIC_COORDS:
        return False, f"No genomic coordinates defined for {case_name}"
    
    coords = CASE_GENOMIC_COORDS[case_name]
    
    try:
        genomic_seq = fetch_genomic_region(
            coords["accession"],
            coords["start"],
            coords["end"],
            coords["strand"],
        )
        
        # Normalize both sequences
        case_norm = case_sequence.upper().replace('T', 'U').replace('\n', '').replace(' ', '')
        genomic_norm = genomic_seq.upper().replace('T', 'U').replace('\n', '').replace(' ', '')
        
        if case_norm == genomic_norm:
            return True, f"Verified: {coords['accession']}:{coords['start']}-{coords['end']}"
        else:
            return False, (
                f"Sequence mismatch for {case_name}\n"
                f"Expected length {len(genomic_norm)}, got {len(case_norm)}\n"
                f"Accession: {coords['accession']}:{coords['start']}-{coords['end']}"
            )
    except Exception as e:
        return False, f"Could not verify {case_name}: {e}"


def fetch_flanking_sequences(
    case_name: str,
    case_sequence: str,
    flank_left: int,
    flank_right: int,
) -> Tuple[str, str]:
    """
    Fetch flanking sequences for window jitter.
    
    Args:
        case_name: Case identifier
        case_sequence: Original case sequence
        flank_left: Number of bases to fetch upstream
        flank_right: Number of bases to fetch downstream
    
    Returns:
        (left_flank, right_flank) tuple
    """
    if case_name not in CASE_GENOMIC_COORDS:
        raise ValueError(f"No coordinates for {case_name}")
    
    coords = CASE_GENOMIC_COORDS[case_name]
    
    # Fetch left flank
    left_start = max(1, coords["start"] - flank_left)
    left_end = coords["start"] - 1
    
    if left_start <= left_end:
        left_flank = fetch_genomic_region(
            coords["accession"],
            left_start,
            left_end,
            coords["strand"],
        )
    else:
        left_flank = ""
    
    # Fetch right flank
    right_start = coords["end"] + 1
    right_end = coords["end"] + flank_right
    
    right_flank = fetch_genomic_region(
        coords["accession"],
        right_start,
        right_end,
        coords["strand"],
    )
    
    return left_flank, right_flank


def compute_tier_from_probs(prob_matrix: np.ndarray, i: int, j: int) -> str:
    """
    Compute tier for a base pair.
    
    Args:
        prob_matrix: Probability matrix
        i, j: 0-based pair indices
        
    Returns:
        "FIRM", "SOFT", "FLOPPY", or "UNPAIRED"
    """
    prob = prob_matrix[i, j]
    
    if prob >= 0.85:
        return "FIRM"
    elif prob >= 0.5:
        return "SOFT"
    elif prob > 0:
        return "FLOPPY"
    else:
        return "UNPAIRED"


def compute_per_nucleotide_tiers(
    sequence: str,
    param_set: str = "Turner2004",
    temperature: float = 37.0,
) -> List[str]:
    """
    Compute per-nucleotide tier classification using MEA structure.
    
    FoldTrust tier logic:
    - Positions paired in the MEA structure get their pair's tier
    - Tier of a pair (i,j) based on p(i,j): FIRM >= 0.85, SOFT >= 0.5, FLOPPY < 0.5
    - Unpaired positions are UNPAIRED
    
    Returns:
        List of tier labels, one per nucleotide
    """
    prob_matrix = compute_pair_probs_with_params(sequence, param_set, temperature)
    n = len(sequence)
    
    # Compute MEA structure
    mea_pairs = compute_mea_structure(prob_matrix)
    
    # Initialize all positions as UNPAIRED
    tiers = ["UNPAIRED"] * n
    
    # Assign tiers based on MEA pairing
    for i, j in mea_pairs:
        pair_prob = prob_matrix[i, j]
        
        if pair_prob >= 0.85:
            tier = "FIRM"
        elif pair_prob >= 0.5:
            tier = "SOFT"
        else:  # pair_prob < 0.5
            tier = "FLOPPY"
        
        # Both positions in the pair get the same tier
        tiers[i] = tier
        tiers[j] = tier
    
    return tiers


def compute_mea_structure(prob_matrix: np.ndarray, gamma: float = 1.0) -> List[Tuple[int, int]]:
    """
    Compute Maximum Expected Accuracy (MEA) structure from probability matrix.
    
    Simple MEA implementation using greedy approach:
    - Score each pair as: 2*gamma*p(i,j) - (sum unpaired penalties)
    - Greedily select highest-scoring compatible pairs
    
    Args:
        prob_matrix: NxN probability matrix
        gamma: Weight for base pairs vs unpaired (default 1.0)
        
    Returns:
        List of (i, j) pairs (0-based)
    """
    n = prob_matrix.shape[0]
    
    # Get all candidate pairs sorted by probability
    candidates = []
    for i in range(n):
        for j in range(i + 1, n):
            if prob_matrix[i, j] > 0:
                candidates.append((prob_matrix[i, j], i, j))
    
    candidates.sort(reverse=True)
    
    # Greedily select compatible pairs
    pairs = []
    paired = set()
    
    for prob, i, j in candidates:
        if i not in paired and j not in paired:
            # Check pseudoknot constraint
            is_compatible = True
            for pi, pj in pairs:
                # Check if (i,j) crosses (pi, pj)
                if (pi < i < pj < j) or (i < pi < j < pj):
                    is_compatible = False
                    break
            
            if is_compatible:
                pairs.append((i, j))
                paired.add(i)
                paired.add(j)
    
    return sorted(pairs)


def basepair_jaccard(pairs1: List[Tuple[int, int]], pairs2: List[Tuple[int, int]]) -> float:
    """
    Compute Jaccard index between two sets of base pairs.
    
    Args:
        pairs1, pairs2: Lists of (i, j) tuples
        
    Returns:
        Jaccard index (intersection / union)
    """
    set1 = set(pairs1)
    set2 = set(pairs2)
    
    if len(set1) == 0 and len(set2) == 0:
        return 1.0
    
    intersection = len(set1 & set2)
    union = len(set1 | set2)
    
    return intersection / union if union > 0 else 0.0


def compute_unpaired_probs(prob_matrix: np.ndarray) -> np.ndarray:
    """
    Compute per-nucleotide unpaired probability.
    
    unpaired[i] = 1 - sum_j p(i,j)
    
    Args:
        prob_matrix: NxN probability matrix
        
    Returns:
        Vector of unpaired probabilities
    """
    paired_probs = np.sum(prob_matrix, axis=1)
    unpaired_probs = 1.0 - paired_probs
    
    # Clamp to [0, 1] (handle numerical precision)
    unpaired_probs = np.clip(unpaired_probs, 0.0, 1.0)
    
    return unpaired_probs


def analyze_robustness_condition(
    sequence: str,
    baseline_tiers: List[str],
    baseline_mea: List[Tuple[int, int]],
    baseline_mfe: List[Tuple[int, int]],
    baseline_unpaired: np.ndarray,
    baseline_firm_pairs: set,
    baseline_floppy_pairs: set,
    param_set: str = "Turner2004",
    temperature: float = 37.0,
) -> Dict:
    """
    Analyze single robustness condition and compute metrics.
    
    Args:
        sequence: RNA sequence
        baseline_*: Baseline reference data
        param_set: Parameter set name
        temperature: Temperature in Celsius
        
    Returns:
        Dictionary of robustness metrics
    """
    # Compute condition structure and probabilities
    try:
        structure, mfe_energy = fold_with_params(sequence, param_set, temperature)
        prob_matrix = compute_pair_probs_with_params(sequence, param_set, temperature)
    except ViennaRNAError as e:
        return {
            "error": str(e),
            "tier_agreement": 0.0,
            "mea_jaccard": 0.0,
            "mfe_jaccard": 0.0,
            "unpaired_spearman": 0.0,
            "firm_retention": 0.0,
            "floppy_retention": 0.0,
            "mfe_energy": None,
        }
    
    # 1. Per-nucleotide tier agreement
    condition_tiers = compute_per_nucleotide_tiers(sequence, param_set, temperature)
    tier_agreement = sum(
        1 for b, c in zip(baseline_tiers, condition_tiers) if b == c
    ) / len(baseline_tiers)
    
    # 2. MEA structure Jaccard
    condition_mea = compute_mea_structure(prob_matrix)
    mea_jaccard = basepair_jaccard(baseline_mea, condition_mea)
    
    # 3. MFE structure Jaccard
    condition_mfe = []
    stack = []
    for i, char in enumerate(structure):
        if char == '(':
            stack.append(i)
        elif char == ')' and stack:
            j = stack.pop()
            condition_mfe.append((j, i))
    mfe_jaccard = basepair_jaccard(baseline_mfe, condition_mfe)
    
    # 4. Unpaired probability Spearman correlation
    condition_unpaired = compute_unpaired_probs(prob_matrix)
    if len(condition_unpaired) == len(baseline_unpaired):
        unpaired_corr, _ = spearmanr(baseline_unpaired, condition_unpaired)
        if np.isnan(unpaired_corr):
            unpaired_corr = 0.0
    else:
        unpaired_corr = 0.0
    
    # 5. FIRM pair retention
    condition_firm_pairs = set()
    for i in range(len(sequence)):
        for j in range(i + 1, len(sequence)):
            if prob_matrix[i, j] >= 0.85:
                condition_firm_pairs.add((i, j))
    
    if len(baseline_firm_pairs) > 0:
        firm_retention = len(baseline_firm_pairs & condition_firm_pairs) / len(baseline_firm_pairs)
    else:
        firm_retention = 1.0  # No FIRM pairs in baseline
    
    # 6. FLOPPY pair retention
    condition_floppy_pairs = set()
    for i in range(len(sequence)):
        for j in range(i + 1, len(sequence)):
            if 0 < prob_matrix[i, j] < 0.5:
                condition_floppy_pairs.add((i, j))
    
    if len(baseline_floppy_pairs) > 0:
        floppy_retention = len(baseline_floppy_pairs & condition_floppy_pairs) / len(baseline_floppy_pairs)
    else:
        floppy_retention = 1.0  # No FLOPPY pairs in baseline
    
    return {
        "tier_agreement": tier_agreement,
        "mea_jaccard": mea_jaccard,
        "mfe_jaccard": mfe_jaccard,
        "unpaired_spearman": unpaired_corr,
        "firm_retention": firm_retention,
        "floppy_retention": floppy_retention,
        "mfe_energy": mfe_energy,
    }


def analyze_case_robustness(
    case_name: str,
    sequence: str,
    max_jitter: int = 25,
) -> Dict:
    """
    Perform comprehensive robustness analysis for one case.
    
    Tests:
    1. Parameter sets: Andronescu2007, Langdon2018 (vs Turner2004 baseline)
    2. Temperatures: 24°C, 30°C, 42°C, 45°C (vs 37°C baseline)
    3. Window jitter: extend/shrink each side by 10 and 25 nt
    
    Args:
        case_name: Case identifier
        sequence: RNA sequence
        max_jitter: Maximum jitter in nt (default 25)
        
    Returns:
        Dictionary with robustness results for all conditions
    """
    print(f"\n=== {case_name} ===")
    print(f"Length: {len(sequence)} nt")
    
    # Compute baseline (Turner2004, 37°C, original window)
    print("Computing baseline...")
    baseline_tiers = compute_per_nucleotide_tiers(sequence, "Turner2004", 37.0)
    baseline_prob_matrix = compute_pair_probs_with_params(sequence, "Turner2004", 37.0)
    baseline_mea = compute_mea_structure(baseline_prob_matrix)
    
    baseline_structure, baseline_mfe_energy = fold_with_params(sequence, "Turner2004", 37.0)
    baseline_mfe_pairs = []
    stack = []
    for i, char in enumerate(baseline_structure):
        if char == '(':
            stack.append(i)
        elif char == ')' and stack:
            j = stack.pop()
            baseline_mfe_pairs.append((j, i))
    
    baseline_unpaired = compute_unpaired_probs(baseline_prob_matrix)
    
    # Extract FIRM and FLOPPY pairs from baseline
    baseline_firm_pairs = set()
    baseline_floppy_pairs = set()
    for i in range(len(sequence)):
        for j in range(i + 1, len(sequence)):
            prob = baseline_prob_matrix[i, j]
            if prob >= 0.85:
                baseline_firm_pairs.add((i, j))
            elif 0 < prob < 0.5:
                baseline_floppy_pairs.add((i, j))
    
    results = {
        "case_name": case_name,
        "length": len(sequence),
        "baseline_mfe_energy": baseline_mfe_energy,
        "baseline_firm_pairs": len(baseline_firm_pairs),
        "baseline_floppy_pairs": len(baseline_floppy_pairs),
        "conditions": {},
    }
    
    # 1. Parameter sets
    print("Testing parameter sets...")
    for param_set in ["Andronescu2007", "Langdon2018"]:
        print(f"  {param_set}...")
        metrics = analyze_robustness_condition(
            sequence,
            baseline_tiers,
            baseline_mea,
            baseline_mfe_pairs,
            baseline_unpaired,
            baseline_firm_pairs,
            baseline_floppy_pairs,
            param_set=param_set,
            temperature=37.0,
        )
        results["conditions"][f"param_{param_set}"] = metrics
    
    # 2. Temperatures
    print("Testing temperatures...")
    for temp in [24.0, 30.0, 42.0, 45.0]:
        print(f"  {temp}°C...")
        metrics = analyze_robustness_condition(
            sequence,
            baseline_tiers,
            baseline_mea,
            baseline_mfe_pairs,
            baseline_unpaired,
            baseline_firm_pairs,
            baseline_floppy_pairs,
            param_set="Turner2004",
            temperature=temp,
        )
        results["conditions"][f"temp_{temp}C"] = metrics
    
    # 3. Window jitter
    print("Fetching flanking sequences for jitter...")
    try:
        left_flank, right_flank = fetch_flanking_sequences(case_name, sequence, max_jitter, max_jitter)
        
        print("Testing window jitter...")
        # Extend left by 10
        jitter_seq = left_flank[-10:] + sequence if len(left_flank) >= 10 else sequence
        metrics = analyze_robustness_condition(
            jitter_seq,
            baseline_tiers,
            baseline_mea,
            baseline_mfe_pairs,
            baseline_unpaired,
            baseline_firm_pairs,
            baseline_floppy_pairs,
            param_set="Turner2004",
            temperature=37.0,
        )
        results["conditions"]["jitter_extend_left_10"] = metrics
        
        # Extend right by 10
        jitter_seq = sequence + right_flank[:10] if len(right_flank) >= 10 else sequence
        metrics = analyze_robustness_condition(
            jitter_seq,
            baseline_tiers,
            baseline_mea,
            baseline_mfe_pairs,
            baseline_unpaired,
            baseline_firm_pairs,
            baseline_floppy_pairs,
            param_set="Turner2004",
            temperature=37.0,
        )
        results["conditions"]["jitter_extend_right_10"] = metrics
        
        # Extend both by 10
        jitter_seq = (left_flank[-10:] if len(left_flank) >= 10 else "") + sequence + (right_flank[:10] if len(right_flank) >= 10 else "")
        metrics = analyze_robustness_condition(
            jitter_seq,
            baseline_tiers,
            baseline_mea,
            baseline_mfe_pairs,
            baseline_unpaired,
            baseline_firm_pairs,
            baseline_floppy_pairs,
            param_set="Turner2004",
            temperature=37.0,
        )
        results["conditions"]["jitter_extend_both_10"] = metrics
        
        # Shrink left by 10
        jitter_seq = sequence[10:] if len(sequence) > 10 else sequence
        metrics = analyze_robustness_condition(
            jitter_seq,
            baseline_tiers[10:] if len(baseline_tiers) > 10 else baseline_tiers,
            baseline_mea,
            baseline_mfe_pairs,
            baseline_unpaired[10:] if len(baseline_unpaired) > 10 else baseline_unpaired,
            baseline_firm_pairs,
            baseline_floppy_pairs,
            param_set="Turner2004",
            temperature=37.0,
        )
        results["conditions"]["jitter_shrink_left_10"] = metrics
        
        # Shrink right by 10
        jitter_seq = sequence[:-10] if len(sequence) > 10 else sequence
        metrics = analyze_robustness_condition(
            jitter_seq,
            baseline_tiers[:-10] if len(baseline_tiers) > 10 else baseline_tiers,
            baseline_mea,
            baseline_mfe_pairs,
            baseline_unpaired[:-10] if len(baseline_unpaired) > 10 else baseline_unpaired,
            baseline_firm_pairs,
            baseline_floppy_pairs,
            param_set="Turner2004",
            temperature=37.0,
        )
        results["conditions"]["jitter_shrink_right_10"] = metrics
        
        # Extend left by 25
        jitter_seq = left_flank[-25:] + sequence if len(left_flank) >= 25 else sequence
        metrics = analyze_robustness_condition(
            jitter_seq,
            baseline_tiers,
            baseline_mea,
            baseline_mfe_pairs,
            baseline_unpaired,
            baseline_firm_pairs,
            baseline_floppy_pairs,
            param_set="Turner2004",
            temperature=37.0,
        )
        results["conditions"]["jitter_extend_left_25"] = metrics
        
        # Extend right by 25
        jitter_seq = sequence + right_flank[:25] if len(right_flank) >= 25 else sequence
        metrics = analyze_robustness_condition(
            jitter_seq,
            baseline_tiers,
            baseline_mea,
            baseline_mfe_pairs,
            baseline_unpaired,
            baseline_firm_pairs,
            baseline_floppy_pairs,
            param_set="Turner2004",
            temperature=37.0,
        )
        results["conditions"]["jitter_extend_right_25"] = metrics
        
        # Extend both by 25
        jitter_seq = (left_flank[-25:] if len(left_flank) >= 25 else "") + sequence + (right_flank[:25] if len(right_flank) >= 25 else "")
        metrics = analyze_robustness_condition(
            jitter_seq,
            baseline_tiers,
            baseline_mea,
            baseline_mfe_pairs,
            baseline_unpaired,
            baseline_firm_pairs,
            baseline_floppy_pairs,
            param_set="Turner2004",
            temperature=37.0,
        )
        results["conditions"]["jitter_extend_both_25"] = metrics
        
        # Shrink left by 25
        jitter_seq = sequence[25:] if len(sequence) > 25 else sequence
        metrics = analyze_robustness_condition(
            jitter_seq,
            baseline_tiers[25:] if len(baseline_tiers) > 25 else baseline_tiers,
            baseline_mea,
            baseline_mfe_pairs,
            baseline_unpaired[25:] if len(baseline_unpaired) > 25 else baseline_unpaired,
            baseline_firm_pairs,
            baseline_floppy_pairs,
            param_set="Turner2004",
            temperature=37.0,
        )
        results["conditions"]["jitter_shrink_left_25"] = metrics
        
        # Shrink right by 25
        jitter_seq = sequence[:-25] if len(sequence) > 25 else sequence
        metrics = analyze_robustness_condition(
            jitter_seq,
            baseline_tiers[:-25] if len(baseline_tiers) > 25 else baseline_tiers,
            baseline_mea,
            baseline_mfe_pairs,
            baseline_unpaired[:-25] if len(baseline_unpaired) > 25 else baseline_unpaired,
            baseline_firm_pairs,
            baseline_floppy_pairs,
            param_set="Turner2004",
            temperature=37.0,
        )
        results["conditions"]["jitter_shrink_right_25"] = metrics
        
    except Exception as e:
        print(f"  Warning: Could not test jitter: {e}")
        results["jitter_error"] = str(e)
    
    return results


def run_robustness_analysis(
    output_dir: Path,
    case_dirs: List[Path],
) -> Dict:
    """
    Run complete robustness analysis for all cases.
    
    Args:
        output_dir: Output directory for results
        case_dirs: List of case directories
        
    Returns:
        Summary dictionary
    """
    print("\n=== Layer 5: Robustness Analysis ===\n")
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Verify ViennaRNA Python is available
    try:
        from foldtrust.vienna import HAS_RNA_PYTHON
        if not HAS_RNA_PYTHON:
            raise ImportError("ViennaRNA Python bindings not available")
    except ImportError as e:
        print(f"ERROR: {e}")
        print("Install with: pip install ViennaRNA")
        return {"error": "ViennaRNA Python bindings required"}
    
    # Process each case
    all_results = []
    
    for case_dir in case_dirs:
        case_name = case_dir.name
        sequence_file = case_dir / "sequence.fa"
        
        if not sequence_file.exists():
            print(f"Skipping {case_name}: no sequence.fa")
            continue
        
        _, sequence = read_fasta(sequence_file)
        
        # Verify sequence against genomic coordinates
        verified, msg = verify_case_sequence(case_name, sequence)
        print(f"\n{case_name} sequence verification: {msg}")
        
        if not verified:
            print(f"WARNING: Sequence verification failed for {case_name}")
            print(f"         {msg}")
            print(f"         Continuing anyway...")
        
        # Run robustness analysis
        try:
            case_results = analyze_case_robustness(case_name, sequence)
            all_results.append(case_results)
        except Exception as e:
            print(f"ERROR analyzing {case_name}: {e}")
            import traceback
            traceback.print_exc()
            continue
    
    # Save detailed results
    results_file = output_dir / "robustness_detailed.json"
    with open(results_file, 'w') as f:
        json.dump(all_results, f, indent=2)
    print(f"\n✓ Saved detailed results to {results_file}")
    
    # Create summary tables
    create_summary_tables(all_results, output_dir)
    
    # Generate figures
    create_robustness_figures(all_results, output_dir)
    
    # Create per-region stability table
    create_region_stability_table(all_results, output_dir)
    
    return {
        "n_cases": len(all_results),
        "output_dir": str(output_dir),
    }


def create_summary_tables(all_results: List[Dict], output_dir: Path) -> None:
    """Create CSV summary tables of robustness metrics."""
    
    # Main robustness table
    rows = []
    for result in all_results:
        case_name = result["case_name"]
        
        for condition_name, metrics in result["conditions"].items():
            if "error" in metrics:
                continue
            
            rows.append({
                "case": case_name,
                "condition": condition_name,
                "tier_agreement": metrics["tier_agreement"],
                "mea_jaccard": metrics["mea_jaccard"],
                "mfe_jaccard": metrics["mfe_jaccard"],
                "unpaired_spearman": metrics["unpaired_spearman"],
                "firm_retention": metrics["firm_retention"],
                "floppy_retention": metrics["floppy_retention"],
                "mfe_energy": metrics.get("mfe_energy"),
            })
    
    df = pd.DataFrame(rows)
    csv_path = output_dir / "robustness_summary.csv"
    df.to_csv(csv_path, index=False, float_format="%.4f")
    print(f"✓ Saved summary table to {csv_path}")
    
    # MFE energy table (parameter sets)
    energy_rows = []
    for result in all_results:
        case_name = result["case_name"]
        energy_row = {
            "case": case_name,
            "Turner2004": result["baseline_mfe_energy"],
        }
        
        for param_set in ["Andronescu2007", "Langdon2018"]:
            cond_name = f"param_{param_set}"
            if cond_name in result["conditions"]:
                energy_row[param_set] = result["conditions"][cond_name].get("mfe_energy")
        
        energy_rows.append(energy_row)
    
    energy_df = pd.DataFrame(energy_rows)
    energy_csv = output_dir / "mfe_energies_by_params.csv"
    energy_df.to_csv(energy_csv, index=False, float_format="%.2f")
    print(f"✓ Saved MFE energies to {energy_csv}")


def create_robustness_figures(all_results: List[Dict], output_dir: Path) -> None:
    """Generate robustness figures."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import seaborn as sns
    
    figures_dir = output_dir.parent / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    
    # Prepare data for heatmaps
    cases = [r["case_name"] for r in all_results]
    
    # Tier agreement heatmap
    tier_data = []
    mea_data = []
    conditions_order = []
    
    # Extract conditions in consistent order
    if all_results:
        first_conditions = list(all_results[0]["conditions"].keys())
        conditions_order = sorted(first_conditions)
    
    for result in all_results:
        tier_row = []
        mea_row = []
        
        for cond in conditions_order:
            if cond in result["conditions"]:
                metrics = result["conditions"][cond]
                tier_row.append(metrics.get("tier_agreement", 0))
                mea_row.append(metrics.get("mea_jaccard", 0))
            else:
                tier_row.append(0)
                mea_row.append(0)
        
        tier_data.append(tier_row)
        mea_data.append(mea_row)
    
    # Shorten condition names for display
    condition_labels = []
    for cond in conditions_order:
        if cond.startswith("param_"):
            label = cond.replace("param_", "")
        elif cond.startswith("temp_"):
            label = cond.replace("temp_", "").replace("C", "°C")
        elif cond.startswith("jitter_"):
            label = cond.replace("jitter_", "").replace("_", " ")
        else:
            label = cond
        condition_labels.append(label)
    
    # Tier agreement heatmap
    fig, ax = plt.subplots(figsize=(12, 6))
    sns.heatmap(
        tier_data,
        xticklabels=condition_labels,
        yticklabels=cases,
        annot=True,
        fmt=".2f",
        cmap="RdYlGn",
        vmin=0.0,
        vmax=1.0,
        ax=ax,
        cbar_kws={"label": "Tier Agreement"},
    )
    ax.set_xlabel("Condition")
    ax.set_ylabel("Case")
    ax.set_title("Per-Nucleotide Tier Agreement Across Conditions")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    tier_fig_path = figures_dir / "layer5_tier_agreement_heatmap.png"
    plt.savefig(tier_fig_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"✓ Saved tier agreement heatmap to {tier_fig_path}")
    
    # MEA Jaccard heatmap
    fig, ax = plt.subplots(figsize=(12, 6))
    sns.heatmap(
        mea_data,
        xticklabels=condition_labels,
        yticklabels=cases,
        annot=True,
        fmt=".2f",
        cmap="RdYlGn",
        vmin=0.0,
        vmax=1.0,
        ax=ax,
        cbar_kws={"label": "MEA Jaccard"},
    )
    ax.set_xlabel("Condition")
    ax.set_ylabel("Case")
    ax.set_title("MEA Structure Jaccard Index Across Conditions")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    mea_fig_path = figures_dir / "layer5_mea_jaccard_heatmap.png"
    plt.savefig(mea_fig_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"✓ Saved MEA Jaccard heatmap to {mea_fig_path}")


def create_region_stability_table(all_results: List[Dict], output_dir: Path) -> None:
    """
    Create per-region stability table for later synthesis.
    
    For each stem/helix in baseline structure, compute:
    - Coordinates
    - Tier
    - Mean pair probability
    - Stability across all conditions (mean tier agreement)
    """
    rows = []
    
    for result in all_results:
        case_name = result["case_name"]
        
        # We don't have direct access to baseline stems here,
        # so we'll create a simplified version showing overall stability
        
        # Compute mean metrics across all conditions
        tier_agreements = []
        for metrics in result["conditions"].values():
            if "error" not in metrics:
                tier_agreements.append(metrics.get("tier_agreement", 0))
        
        if tier_agreements:
            mean_stability = np.mean(tier_agreements)
            std_stability = np.std(tier_agreements)
        else:
            mean_stability = 0.0
            std_stability = 0.0
        
        rows.append({
            "case": case_name,
            "length": result["length"],
            "baseline_firm_pairs": result["baseline_firm_pairs"],
            "baseline_floppy_pairs": result["baseline_floppy_pairs"],
            "mean_tier_agreement": mean_stability,
            "std_tier_agreement": std_stability,
        })
    
    df = pd.DataFrame(rows)
    csv_path = output_dir / "region_stability.csv"
    df.to_csv(csv_path, index=False, float_format="%.4f")
    print(f"✓ Saved region stability table to {csv_path}")
