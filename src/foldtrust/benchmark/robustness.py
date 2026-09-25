"""Robustness analysis: parameter sets, temperature, and window boundaries."""

import json
from pathlib import Path
from typing import Dict, List, Tuple
import pandas as pd
import numpy as np

try:
    import RNA
    HAS_RNA_PYTHON = True
except ImportError:
    HAS_RNA_PYTHON = False

from foldtrust.vienna import compute_pair_probabilities, parse_stems
from foldtrust.utils import read_fasta


def fold_with_params(sequence: str, param_set: str = "Turner2004", temperature: float = 37.0) -> Tuple[str, float, np.ndarray]:
    """
    Fold RNA with specific parameter set and temperature using ViennaRNA Python.
    
    Args:
        sequence: RNA sequence
        param_set: One of "Turner2004", "Turner1999", "Andronescu2007", "Langdon2018"
        temperature: Temperature in Celsius
        
    Returns:
        (structure, mfe_energy, pair_probability_matrix)
    """
    if not HAS_RNA_PYTHON:
        raise RuntimeError("ViennaRNA Python module required. Install: pip install ViennaRNA")
    
    # Load parameter set
    if param_set == "Turner2004":
        RNA.params_load_RNA_Turner2004()
    elif param_set == "Turner1999":
        RNA.params_load_RNA_Turner1999()
    elif param_set == "Andronescu2007":
        RNA.params_load_RNA_Andronescu2007()
    elif param_set == "Langdon2018":
        RNA.params_load_RNA_Langdon2018()
    else:
        raise ValueError(f"Unknown parameter set: {param_set}")
    
    # Create fold compound with temperature
    md = RNA.md()
    md.temperature = temperature
    fc = RNA.fold_compound(sequence, md)
    
    # Compute MFE
    structure, mfe = fc.mfe()
    
    # Compute partition function and get base pair probabilities
    fc.pf()
    bpp = fc.bpp()
    
    # Convert to numpy matrix (bpp is 1-indexed, returns n+1 x n+1)
    n = len(sequence)
    prob_matrix = np.zeros((n, n))
    
    for i in range(n):
        for j in range(n):
            if i < len(bpp) and j < len(bpp[i]):
                prob_matrix[i, j] = bpp[i][j]
    
    return structure, mfe, prob_matrix


def compare_tier_classifications(stems1: List[Dict], stems2: List[Dict]) -> Dict:
    """
    Compare tier classifications between two stem sets.
    
    Returns:
        Dictionary with stability metrics
    """
    # Build mapping of pair -> tier
    tier_map1 = {}
    for stem in stems1:
        for pair in stem["pairs"]:
            tier_map1[pair] = stem["flag"]
    
    tier_map2 = {}
    for stem in stems2:
        for pair in stem["pairs"]:
            tier_map2[pair] = stem["flag"]
    
    # Find common pairs
    common_pairs = set(tier_map1.keys()) & set(tier_map2.keys())
    
    if len(common_pairs) == 0:
        return {
            "common_pairs": 0,
            "changed_pairs": 0,
            "stability": 0.0,
            "changes": {}
        }
    
    # Count changes
    changes = {"firm_to_soft": 0, "firm_to_floppy": 0, "soft_to_firm": 0, 
               "soft_to_floppy": 0, "floppy_to_firm": 0, "floppy_to_soft": 0}
    
    changed = 0
    for pair in common_pairs:
        t1, t2 = tier_map1[pair], tier_map2[pair]
        if t1 != t2:
            changed += 1
            key = f"{t1}_to_{t2}"
            if key in changes:
                changes[key] += 1
    
    return {
        "common_pairs": len(common_pairs),
        "changed_pairs": changed,
        "stability": 1.0 - (changed / len(common_pairs)) if common_pairs else 0.0,
        "changes": changes
    }


def run_parameter_set_comparison(
    output_dir: Path,
    case_dirs: List[Path],
) -> Dict:
    """
    Compare tier classifications across different parameter sets.
    
    Args:
        output_dir: Output directory
        case_dirs: List of case directories
        
    Returns:
        Dictionary with comparison results
    """
    if not HAS_RNA_PYTHON:
        print("  ViennaRNA Python module not available. Skipping parameter set comparison.")
        return {}
    
    print("\n=== Parameter Set Comparison ===\n")
    
    param_sets = ["Turner2004", "Andronescu2007", "Langdon2018"]
    results = []
    
    for case_dir in case_dirs:
        case_name = case_dir.name
        sequence_file = case_dir / "sequence.fa"
        
        if not sequence_file.exists():
            continue
        
        print(f"Analyzing {case_name}...")
        
        _, sequence = read_fasta(sequence_file)
        
        # Fold with each parameter set
        param_stems = {}
        for param_set in param_sets:
            try:
                structure, mfe, prob_matrix = fold_with_params(sequence, param_set, 37.0)
                stems = parse_stems(structure, prob_matrix)
                param_stems[param_set] = stems
            except Exception as e:
                print(f"  Error with {param_set}: {e}")
                continue
        
        # Compare Turner2004 (baseline) to others
        if "Turner2004" in param_stems:
            baseline = param_stems["Turner2004"]
            
            for param_set in ["Andronescu2007", "Langdon2018"]:
                if param_set in param_stems:
                    comparison = compare_tier_classifications(baseline, param_stems[param_set])
                    
                    results.append({
                        "case": case_name,
                        "length": len(sequence),
                        "comparison": f"Turner2004_vs_{param_set}",
                        "common_pairs": comparison["common_pairs"],
                        "changed_pairs": comparison["changed_pairs"],
                        "stability": comparison["stability"],
                        **comparison["changes"]
                    })
    
    df = pd.DataFrame(results)
    
    csv_path = output_dir / "parameter_set_comparison.csv"
    df.to_csv(csv_path, index=False)
    print(f"\n✓ Saved to {csv_path}")
    
    summary = {
        "n_cases": len(set(df["case"])) if len(df) > 0 else 0,
        "mean_stability": float(df["stability"].mean()) if len(df) > 0 else 0.0,
        "parameter_sets": param_sets,
    }
    
    return summary


def run_temperature_sweep(
    output_dir: Path,
    case_dirs: List[Path],
    temperatures: List[float] = [24.0, 37.0, 42.0],
) -> Dict:
    """
    Analyze tier stability across temperatures.
    
    Args:
        output_dir: Output directory
        case_dirs: List of case directories
        temperatures: List of temperatures in Celsius
        
    Returns:
        Dictionary with temperature sweep results
    """
    if not HAS_RNA_PYTHON:
        print("  ViennaRNA Python module not available. Skipping temperature sweep.")
        return {}
    
    print("\n=== Temperature Sweep ===\n")
    
    results = []
    
    for case_dir in case_dirs:
        case_name = case_dir.name
        sequence_file = case_dir / "sequence.fa"
        
        if not sequence_file.exists():
            continue
        
        print(f"Analyzing {case_name}...")
        
        _, sequence = read_fasta(sequence_file)
        
        # Fold at each temperature
        temp_stems = {}
        for temp in temperatures:
            try:
                structure, mfe, prob_matrix = fold_with_params(sequence, "Turner2004", temp)
                stems = parse_stems(structure, prob_matrix)
                temp_stems[temp] = stems
                print(f"  {temp}°C: {len(stems)} stems")
            except Exception as e:
                print(f"  Error at {temp}°C: {e}")
                continue
        
        # Compare 37°C (baseline) to other temperatures
        if 37.0 in temp_stems:
            baseline = temp_stems[37.0]
            
            for temp in temperatures:
                if temp != 37.0 and temp in temp_stems:
                    comparison = compare_tier_classifications(baseline, temp_stems[temp])
                    
                    results.append({
                        "case": case_name,
                        "length": len(sequence),
                        "baseline_temp": 37.0,
                        "comparison_temp": temp,
                        "common_pairs": comparison["common_pairs"],
                        "changed_pairs": comparison["changed_pairs"],
                        "stability": comparison["stability"],
                        **comparison["changes"]
                    })
    
    df = pd.DataFrame(results)
    
    csv_path = output_dir / "temperature_sweep.csv"
    df.to_csv(csv_path, index=False)
    print(f"\n✓ Saved to {csv_path}")
    
    summary = {
        "n_cases": len(set(df["case"])) if len(df) > 0 else 0,
        "temperatures": temperatures,
        "mean_stability_24C": float(df[df["comparison_temp"] == 24.0]["stability"].mean()) if len(df[df["comparison_temp"] == 24.0]) > 0 else 0.0,
        "mean_stability_42C": float(df[df["comparison_temp"] == 42.0]["stability"].mean()) if len(df[df["comparison_temp"] == 42.0]) > 0 else 0.0,
    }
    
    return summary


def run_window_jitter_analysis(
    output_dir: Path,
    case_data: List[Dict],
) -> Dict:
    """
    Analyze tier stability under window boundary perturbations.
    
    Args:
        output_dir: Output directory
        case_data: List of dicts with 'name', 'sequence', 'extended_sequence'
        
    Returns:
        Dictionary with jitter analysis results
    """
    if not HAS_RNA_PYTHON:
        print("  ViennaRNA Python module not available. Skipping window jitter.")
        return {}
    
    print("\n=== Window Jitter Analysis ===\n")
    
    jitters = [-25, -10, 10, 25]
    results = []
    
    for case in case_data:
        case_name = case["name"]
        baseline_seq = case["sequence"]
        extended_seq = case.get("extended_sequence")
        
        if not extended_seq:
            print(f"Skipping {case_name}: no extended sequence")
            continue
        
        print(f"Analyzing {case_name}...")
        
        # Find baseline window in extended sequence
        baseline_start = extended_seq.find(baseline_seq)
        if baseline_start < 0:
            print(f"  Error: baseline sequence not found in extended sequence")
            continue
        
        # Fold baseline
        try:
            base_struct, base_mfe, base_prob = fold_with_params(baseline_seq, "Turner2004", 37.0)
            base_stems = parse_stems(base_struct, base_prob)
        except Exception as e:
            print(f"  Error folding baseline: {e}")
            continue
        
        # Try each jitter
        for jitter in jitters:
            new_start = baseline_start + jitter
            new_end = new_start + len(baseline_seq)
            
            if new_start < 0 or new_end > len(extended_seq):
                print(f"  Skipping jitter {jitter:+d}: out of bounds")
                continue
            
            jitter_seq = extended_seq[new_start:new_end]
            
            try:
                jitter_struct, jitter_mfe, jitter_prob = fold_with_params(jitter_seq, "Turner2004", 37.0)
                jitter_stems = parse_stems(jitter_struct, jitter_prob)
                
                # Note: comparing different sequences, so we compare tier distributions not individual pairs
                base_tiers = {"firm": 0, "soft": 0, "floppy": 0}
                jitter_tiers = {"firm": 0, "soft": 0, "floppy": 0}
                
                for stem in base_stems:
                    base_tiers[stem["flag"]] += 1
                for stem in jitter_stems:
                    jitter_tiers[stem["flag"]] += 1
                
                results.append({
                    "case": case_name,
                    "jitter": jitter,
                    "base_firm": base_tiers["firm"],
                    "base_soft": base_tiers["soft"],
                    "base_floppy": base_tiers["floppy"],
                    "jitter_firm": jitter_tiers["firm"],
                    "jitter_soft": jitter_tiers["soft"],
                    "jitter_floppy": jitter_tiers["floppy"],
                })
                
                print(f"  Jitter {jitter:+3d}nt: {jitter_tiers['firm']}/{jitter_tiers['soft']}/{jitter_tiers['floppy']} (F/S/FL)")
                
            except Exception as e:
                print(f"  Error with jitter {jitter}: {e}")
                continue
    
    df = pd.DataFrame(results)
    
    if len(df) > 0:
        csv_path = output_dir / "window_jitter.csv"
        df.to_csv(csv_path, index=False)
        print(f"\n✓ Saved to {csv_path}")
    
    summary = {
        "n_cases": len(set(df["case"])) if len(df) > 0 else 0,
        "jitters_tested": jitters if len(df) > 0 else [],
    }
    
    return summary


def run_robustness_analysis(
    output_dir: Path,
    case_dirs: List[Path],
    case_extended_sequences: Dict[str, str] = None,
) -> Dict:
    """
    Run complete robustness analysis.
    
    Args:
        output_dir: Output directory
        case_dirs: List of case directories
        case_extended_sequences: Optional dict mapping case names to extended sequences for jitter
        
    Returns:
        Combined summary
    """
    print("\n=== Robustness Analysis ===")
    
    summary = {}
    
    # Parameter set comparison
    param_summary = run_parameter_set_comparison(output_dir, case_dirs)
    summary["parameter_sets"] = param_summary
    
    # Temperature sweep
    temp_summary = run_temperature_sweep(output_dir, case_dirs)
    summary["temperature"] = temp_summary
    
    # Window jitter (if extended sequences provided)
    if case_extended_sequences:
        case_data = []
        for case_dir in case_dirs:
            case_name = case_dir.name
            if case_name in case_extended_sequences:
                _, sequence = read_fasta(case_dir / "sequence.fa")
                case_data.append({
                    "name": case_name,
                    "sequence": sequence,
                    "extended_sequence": case_extended_sequences[case_name]
                })
        
        jitter_summary = run_window_jitter_analysis(output_dir, case_data)
        summary["window_jitter"] = jitter_summary
    
    # Save combined summary
    summary_path = output_dir / "robustness_summary.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)
    
    return summary
