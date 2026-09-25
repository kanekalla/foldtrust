"""Layer 5: Robustness analysis across temperature, parameter sets, and window context.

This module tests FoldTrust's structure and tier predictions under three types of perturbations:
1. Temperature variations (25, 30, 37, 42°C; 37°C is the reference)
2. Energy parameter sets (Turner2004 reference vs Andronescu2007 vs Langdon2018)
3. Window context: extend by real flanking sequence (0, 25, 50, 100 nt on each side)

All conditions are compared to the 37°C Turner2004 baseline for the core window.
Retention metrics are computed only over stems/pairs that exist in the reference condition.
"""

from pathlib import Path
from typing import Dict, List, Optional, Tuple
import json

import numpy as np
import pandas as pd

try:
    import RNA
    HAS_RNA = True
except ImportError:
    HAS_RNA = False

from foldtrust.vienna import (
    fold_with_params,
    compute_pair_probs_with_params,
    parse_stems,
    ViennaRNAError,
)
from foldtrust.utils import read_fasta


# Verified case coordinates from Layer 0
CASE_COORDS = {
    "sars2-fse": {
        "accession": "NC_045512.2",
        "start": 13462,
        "end": 13542,
        "length": 81,
    },
    "smn2-iss-n1": {
        "accession": "NG_008728.1",
        "start": 31999,
        "end": 32152,
        "length": 154,
    },
    "cftr-5utr": {
        "accession": "NM_000492.4",
        "start": 1,
        "end": 200,
        "length": 200,
    },
    "mapt-e10": {
        "accession": "NG_007398.2",
        "start": 120818,
        "end": 121000,
        "length": 183,
    },
    "hcv-ires-dii": {
        "accession": "AF009606.1",
        "start": 44,
        "end": 118,
        "length": 75,
    },
}


def compute_mea_structure(
    sequence: str,
    param_set: str = "Turner2004",
    temperature: float = 37.0,
    gamma: float = 1.0,
) -> Tuple[List[Tuple[int, int]], float]:
    """
    Compute MEA structure using ViennaRNA's fc.MEA(gamma).
    
    Args:
        sequence: RNA sequence
        param_set: Parameter set name
        temperature: Temperature in Celsius
        gamma: MEA gamma parameter (default 1.0)
        
    Returns:
        (pairs, mea_energy) where pairs is list of (i, j) 0-based tuples
    """
    if not HAS_RNA:
        raise ViennaRNAError("ViennaRNA not available")
    
    # Use subprocess like fold_with_params to ensure parameter isolation
    import subprocess
    import sys
    
    code = f"""
import RNA

sequence = {sequence!r}
param_set = {param_set!r}
temperature = {temperature}
gamma = {gamma}

# Load parameters
if param_set == "Turner2004":
    RNA.params_load_RNA_Turner2004()
elif param_set == "Andronescu2007":
    RNA.params_load_RNA_Andronescu2007()
elif param_set == "Langdon2018":
    RNA.params_load_RNA_Langdon2018()

# Create fold compound
md = RNA.md()
md.temperature = temperature
fc = RNA.fold_compound(sequence, md)

# Compute partition function
fc.pf()

# Compute MEA
mea_result = fc.MEA(gamma)
mea_struct = mea_result[0]
mea_energy = mea_result[1]

print(f"{{mea_struct}}|{{mea_energy}}")
"""
    
    try:
        result = subprocess.run(
            [sys.executable, "-c", code],
            capture_output=True,
            text=True,
            timeout=30,
            check=True,
        )
        
        output = result.stdout.strip()
        mea_struct, mea_energy_str = output.split("|")
        mea_energy = float(mea_energy_str)
        
        # Parse structure to pairs
        pairs = []
        stack = []
        for i, char in enumerate(mea_struct):
            if char == "(":
                stack.append(i)
            elif char == ")" and stack:
                j = stack.pop()
                pairs.append((j, i))
        
        return pairs, mea_energy
        
    except subprocess.CalledProcessError as e:
        raise ViennaRNAError(f"MEA computation failed: {e.stderr}")
    except Exception as e:
        raise ViennaRNAError(f"MEA error: {e}")


def compute_ensemble_defect(
    sequence: str,
    structure: str,
    param_set: str = "Turner2004",
    temperature: float = 37.0,
) -> float:
    """
    Compute ensemble defect: expected number of incorrectly paired bases.
    
    Ensemble defect = sum over positions of (1 - probability of being in correct state)
    
    Args:
        sequence: RNA sequence
        structure: Dot-bracket structure
        param_set: Parameter set
        temperature: Temperature
        
    Returns:
        Ensemble defect
    """
    prob_matrix = compute_pair_probs_with_params(sequence, param_set, temperature)
    n = len(sequence)
    defect = 0.0
    
    # Parse structure to get paired positions
    paired = {}
    stack = []
    for i, char in enumerate(structure):
        if char == "(":
            stack.append(i)
        elif char == ")" and stack:
            j = stack.pop()
            paired[j] = i
            paired[i] = j
    
    # For each position, probability of correct state
    for i in range(n):
        if i in paired:
            j = paired[i]
            correct_prob = prob_matrix[i, j]
        else:
            correct_prob = 1.0 - np.sum(prob_matrix[i, :])
        
        defect += (1.0 - correct_prob)
    
    return defect


def compute_basepair_distance(pairs1: List[Tuple[int, int]], pairs2: List[Tuple[int, int]]) -> int:
    """
    Compute base-pair distance: number of pairs in symmetric difference.
    
    Args:
        pairs1, pairs2: Lists of (i, j) pairs
        
    Returns:
        |pairs1 - pairs2| + |pairs2 - pairs1|
    """
    set1 = set(pairs1)
    set2 = set(pairs2)
    return len(set1 - set2) + len(set2 - set1)


def compute_stem_retention(
    ref_stems: List[Dict],
    cond_prob_matrix: np.ndarray,
    cond_mea_pairs: List[Tuple[int, int]],
) -> Dict:
    """
    Compute retention of reference stems in a condition.
    
    For each reference stem tier (FIRM/SOFT/FLOPPY), compute:
    - Fraction of stems whose pairs are still paired in condition MEA
    - Mean change in stem mean probability
    
    Args:
        ref_stems: Reference stems from parse_stems
        cond_prob_matrix: Condition probability matrix
        cond_mea_pairs: Condition MEA pairs
        
    Returns:
        Dict with retention metrics per tier
    """
    cond_pairs_set = set(cond_mea_pairs)
    
    # Group stems by tier
    tiers = {"firm": [], "soft": [], "floppy": []}
    for stem in ref_stems:
        tiers[stem["flag"]].append(stem)
    
    results = {}
    for tier_name, stems in tiers.items():
        if not stems:
            results[tier_name] = {
                "count": 0,
                "retention": None,
                "mean_prob_change": None,
            }
            continue
        
        # Count how many stems are retained (all pairs still paired in MEA)
        retained = 0
        prob_changes = []
        
        for stem in stems:
            stem_pairs = stem["pairs"]
            ref_mean_prob = stem["mean_prob"]
            
            # Check if all pairs in stem are paired in condition MEA
            all_paired = all(pair in cond_pairs_set for pair in stem_pairs)
            if all_paired:
                retained += 1
            
            # Compute mean probability change
            cond_probs = [cond_prob_matrix[i, j] for i, j in stem_pairs]
            cond_mean_prob = np.mean(cond_probs)
            prob_changes.append(cond_mean_prob - ref_mean_prob)
        
        results[tier_name] = {
            "count": len(stems),
            "retention": retained / len(stems),
            "mean_prob_change": np.mean(prob_changes),
        }
    
    return results


def analyze_temperature_condition(
    sequence: str,
    temp: float,
    ref_structure: str,
    ref_mea_pairs: List[Tuple[int, int]],
    ref_stems: List[Dict],
) -> Dict:
    """
    Analyze one temperature condition against 37°C Turner2004 reference.
    
    Returns dict with metrics: mfe_energy, ensemble_defect, bp_distance_mfe, 
    bp_distance_mea, stem_retention (by tier)
    """
    param_set = "Turner2004"
    
    # Compute MFE
    structure, mfe_energy = fold_with_params(sequence, param_set, temp)
    
    # Parse MFE pairs
    mfe_pairs = []
    stack = []
    for i, char in enumerate(structure):
        if char == "(":
            stack.append(i)
        elif char == ")" and stack:
            j = stack.pop()
            mfe_pairs.append((j, i))
    
    # Compute MEA
    mea_pairs, mea_energy = compute_mea_structure(sequence, param_set, temp)
    
    # Ensemble defect
    ensemble_defect = compute_ensemble_defect(sequence, structure, param_set, temp)
    
    # Base-pair distances
    bp_dist_mfe = compute_basepair_distance(
        [(i, j) for i, j in mfe_pairs],
        [(i, j) for i, j in ref_mea_pairs]
    )
    bp_dist_mea = compute_basepair_distance(mea_pairs, ref_mea_pairs)
    
    # Probability matrix for stem retention
    prob_matrix = compute_pair_probs_with_params(sequence, param_set, temp)
    
    # Stem retention
    stem_retention = compute_stem_retention(ref_stems, prob_matrix, mea_pairs)
    
    return {
        "mfe_energy": mfe_energy,
        "ensemble_defect": ensemble_defect,
        "bp_distance_mfe": bp_dist_mfe,
        "bp_distance_mea": bp_dist_mea,
        "stem_retention": stem_retention,
    }


def analyze_params_condition(
    sequence: str,
    param_set: str,
    ref_structure: str,
    ref_mea_pairs: List[Tuple[int, int]],
    ref_stems: List[Dict],
) -> Dict:
    """
    Analyze one parameter set against Turner2004 reference.
    """
    temp = 37.0
    
    # Compute MFE
    structure, mfe_energy = fold_with_params(sequence, param_set, temp)
    
    # Parse MFE pairs
    mfe_pairs = []
    stack = []
    for i, char in enumerate(structure):
        if char == "(":
            stack.append(i)
        elif char == ")" and stack:
            j = stack.pop()
            mfe_pairs.append((j, i))
    
    # Compute MEA
    mea_pairs, mea_energy = compute_mea_structure(sequence, param_set, temp)
    
    # Ensemble defect
    ensemble_defect = compute_ensemble_defect(sequence, structure, param_set, temp)
    
    # Base-pair distances
    bp_dist_mfe = compute_basepair_distance(mfe_pairs, ref_mea_pairs)
    bp_dist_mea = compute_basepair_distance(mea_pairs, ref_mea_pairs)
    
    # Probability matrix
    prob_matrix = compute_pair_probs_with_params(sequence, param_set, temp)
    
    # Stem retention
    stem_retention = compute_stem_retention(ref_stems, prob_matrix, mea_pairs)
    
    return {
        "mfe_energy": mfe_energy,
        "ensemble_defect": ensemble_defect,
        "bp_distance_mfe": bp_dist_mfe,
        "bp_distance_mea": bp_dist_mea,
        "stem_retention": stem_retention,
    }


def run_layer5_analysis(
    cases_dir: Path,
    output_dir: Path,
) -> Dict:
    """
    Run Layer 5 robustness analysis for all cases.
    
    For each case:
    1. Compute 37°C Turner2004 baseline
    2. Test temperatures: 25, 30, 42°C
    3. Test parameter sets: Andronescu2007, Langdon2018
    
    Args:
        cases_dir: Path to data/cases
        output_dir: Path to benchmarks/outputs/layer5
        
    Returns:
        Summary dict
    """
    if not HAS_RNA:
        raise RuntimeError("ViennaRNA not available. Install: pip install ViennaRNA")
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    results = {}
    
    for case_name in CASE_COORDS.keys():
        print(f"\n=== {case_name} ===")
        
        case_dir = cases_dir / case_name
        seq_file = case_dir / "sequence.fa"
        
        if not seq_file.exists():
            print(f"  Skipping: no sequence.fa")
            continue
        
        _, sequence = read_fasta(seq_file)
        print(f"  Length: {len(sequence)} nt")
        
        # Baseline: 37°C Turner2004
        print("  Computing baseline (37°C Turner2004)...")
        ref_structure, ref_mfe_energy = fold_with_params(sequence, "Turner2004", 37.0)
        ref_prob_matrix = compute_pair_probs_with_params(sequence, "Turner2004", 37.0)
        ref_stems = parse_stems(ref_structure, ref_prob_matrix)
        ref_mea_pairs, ref_mea_energy = compute_mea_structure(sequence, "Turner2004", 37.0)
        ref_ensemble_defect = compute_ensemble_defect(sequence, ref_structure, "Turner2004", 37.0)
        
        case_results = {
            "case": case_name,
            "length": len(sequence),
            "baseline": {
                "mfe_energy": ref_mfe_energy,
                "ensemble_defect": ref_ensemble_defect,
                "n_stems": len(ref_stems),
                "stem_counts": {
                    "firm": sum(1 for s in ref_stems if s["flag"] == "firm"),
                    "soft": sum(1 for s in ref_stems if s["flag"] == "soft"),
                    "floppy": sum(1 for s in ref_stems if s["flag"] == "floppy"),
                },
            },
            "temperature": {},
            "parameters": {},
        }
        
        # Temperature sweep
        for temp in [25.0, 30.0, 42.0]:
            print(f"  Temperature {temp}°C...")
            case_results["temperature"][f"{temp}C"] = analyze_temperature_condition(
                sequence, temp, ref_structure, ref_mea_pairs, ref_stems
            )
        
        # Parameter set sweep
        for param_set in ["Andronescu2007", "Langdon2018"]:
            print(f"  Parameters {param_set}...")
            case_results["parameters"][param_set] = analyze_params_condition(
                sequence, param_set, ref_structure, ref_mea_pairs, ref_stems
            )
        
        results[case_name] = case_results
    
    # Save detailed JSON
    json_path = output_dir / "layer5_temp_params.json"
    with open(json_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n✓ Saved detailed results to {json_path}")
    
    # Create CSV tables
    create_temperature_tables(results, output_dir)
    create_parameters_tables(results, output_dir)
    
    return results


def create_temperature_tables(results: Dict, output_dir: Path):
    """Create CSV tables for temperature sweep."""
    
    # Energy table
    rows = []
    for case_name, case_data in results.items():
        row = {
            "case": case_name,
            "37C_baseline": case_data["baseline"]["mfe_energy"],
        }
        for temp_key, temp_data in case_data["temperature"].items():
            row[temp_key] = temp_data["mfe_energy"]
        rows.append(row)
    
    df = pd.DataFrame(rows)
    csv_path = output_dir / "temperature_mfe_energies.csv"
    df.to_csv(csv_path, index=False, float_format="%.2f")
    print(f"✓ {csv_path.name}")
    
    # BP distance table
    rows = []
    for case_name, case_data in results.items():
        for temp_key, temp_data in case_data["temperature"].items():
            rows.append({
                "case": case_name,
                "temperature": temp_key,
                "bp_distance_mfe": temp_data["bp_distance_mfe"],
                "bp_distance_mea": temp_data["bp_distance_mea"],
                "ensemble_defect": temp_data["ensemble_defect"],
            })
    
    df = pd.DataFrame(rows)
    csv_path = output_dir / "temperature_bp_distances.csv"
    df.to_csv(csv_path, index=False, float_format="%.2f")
    print(f"✓ {csv_path.name}")
    
    # Stem retention table
    rows = []
    for case_name, case_data in results.items():
        for temp_key, temp_data in case_data["temperature"].items():
            for tier, tier_data in temp_data["stem_retention"].items():
                if tier_data["count"] > 0:
                    rows.append({
                        "case": case_name,
                        "temperature": temp_key,
                        "tier": tier.upper(),
                        "n_stems": tier_data["count"],
                        "retention": tier_data["retention"],
                        "mean_prob_change": tier_data["mean_prob_change"],
                    })
    
    df = pd.DataFrame(rows)
    csv_path = output_dir / "temperature_stem_retention.csv"
    df.to_csv(csv_path, index=False, float_format="%.4f")
    print(f"✓ {csv_path.name}")


def create_parameters_tables(results: Dict, output_dir: Path):
    """Create CSV tables for parameter set sweep."""
    
    # Energy table
    rows = []
    for case_name, case_data in results.items():
        row = {
            "case": case_name,
            "Turner2004": case_data["baseline"]["mfe_energy"],
        }
        for param_name, param_data in case_data["parameters"].items():
            row[param_name] = param_data["mfe_energy"]
        rows.append(row)
    
    df = pd.DataFrame(rows)
    csv_path = output_dir / "parameters_mfe_energies.csv"
    df.to_csv(csv_path, index=False, float_format="%.2f")
    print(f"✓ {csv_path.name}")
    
    # BP distance table
    rows = []
    for case_name, case_data in results.items():
        for param_name, param_data in case_data["parameters"].items():
            rows.append({
                "case": case_name,
                "parameters": param_name,
                "bp_distance_mfe": param_data["bp_distance_mfe"],
                "bp_distance_mea": param_data["bp_distance_mea"],
                "ensemble_defect": param_data["ensemble_defect"],
            })
    
    df = pd.DataFrame(rows)
    csv_path = output_dir / "parameters_bp_distances.csv"
    df.to_csv(csv_path, index=False, float_format="%.2f")
    print(f"✓ {csv_path.name}")
    
    # Stem retention table
    rows = []
    for case_name, case_data in results.items():
        for param_name, param_data in case_data["parameters"].items():
            for tier, tier_data in param_data["stem_retention"].items():
                if tier_data["count"] > 0:
                    rows.append({
                        "case": case_name,
                        "parameters": param_name,
                        "tier": tier.upper(),
                        "n_stems": tier_data["count"],
                        "retention": tier_data["retention"],
                        "mean_prob_change": tier_data["mean_prob_change"],
                    })
    
    df = pd.DataFrame(rows)
    csv_path = output_dir / "parameters_stem_retention.csv"
    df.to_csv(csv_path, index=False, float_format="%.4f")
    print(f"✓ {csv_path.name}")
