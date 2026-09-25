"""Robustness analysis: parameter sets, temperature, and window boundaries."""

import json
from pathlib import Path
from typing import Dict, List
import pandas as pd
import numpy as np

from foldtrust.vienna import fold_mfe, compute_pair_probabilities, parse_stems
from foldtrust.utils import read_fasta


def run_robustness_analysis(
    output_dir: Path,
    case_dirs: List[Path],
) -> Dict:
    """
    Analyze robustness of reliability tiers under perturbations.
    
    Args:
        output_dir: Directory for outputs
        case_dirs: List of case directories with sequence.fa files
        
    Returns:
        Dictionary with robustness results
    """
    print("\n=== Robustness Analysis ===\n")
    
    results = []
    
    for case_dir in case_dirs:
        case_name = case_dir.name
        sequence_file = case_dir / "sequence.fa"
        
        if not sequence_file.exists():
            continue
        
        print(f"Analyzing {case_name}...")
        
        _, sequence = read_fasta(sequence_file)
        
        baseline_result = analyze_baseline(sequence, case_name)
        
        temp_result = analyze_temperature_robustness(sequence, case_name)
        
        results.append({
            "case_name": case_name,
            "length": len(sequence),
            **baseline_result,
            **temp_result,
        })
    
    df = pd.DataFrame(results)
    
    csv_path = output_dir / "robustness_results.csv"
    df.to_csv(csv_path, index=False)
    
    summary = compute_robustness_summary(df)
    
    summary_path = output_dir / "robustness_summary.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)
    
    generate_robustness_report(df, summary, output_dir)
    
    print(f"\n✓ Completed robustness analysis for {len(results)} cases")
    
    return summary


def analyze_baseline(sequence: str, case_name: str) -> Dict:
    """Analyze baseline tier distribution."""
    
    structure, mfe_energy = fold_mfe(sequence)
    prob_matrix = compute_pair_probabilities(sequence)
    stems = parse_stems(structure, prob_matrix)
    
    firm = sum(1 for s in stems if s["flag"] == "firm")
    soft = sum(1 for s in stems if s["flag"] == "soft")
    floppy = sum(1 for s in stems if s["flag"] == "floppy")
    
    return {
        "baseline_firm": firm,
        "baseline_soft": soft,
        "baseline_floppy": floppy,
        "baseline_total": len(stems),
    }


def analyze_temperature_robustness(sequence: str, case_name: str) -> Dict:
    """
    Analyze tier stability across temperatures.
    
    Note: ViennaRNA CLI doesn't easily support temperature changes without
    recompiling or using the Python bindings. For MVP, we document the 
    limitation and use baseline only.
    """
    
    return {
        "temp_note": "Temperature sweep requires ViennaRNA Python bindings (future work)",
    }


def compute_robustness_summary(df: pd.DataFrame) -> Dict:
    """Compute summary statistics for robustness."""
    
    summary = {
        "n_cases": len(df),
        "baseline_distribution": {
            "firm_mean": float(df["baseline_firm"].mean()),
            "soft_mean": float(df["baseline_soft"].mean()),
            "floppy_mean": float(df["baseline_floppy"].mean()),
            "total_mean": float(df["baseline_total"].mean()),
        },
        "note": "Temperature and parameter set sweep requires ViennaRNA Python library (future enhancement)",
    }
    
    return summary


def generate_robustness_report(df: pd.DataFrame, summary: Dict, output_dir: Path) -> None:
    """Generate markdown report for robustness analysis."""
    
    report_lines = [
        "# Robustness Analysis",
        "",
        "## Baseline Tier Distribution",
        "",
        f"- **N cases:** {summary['n_cases']}",
        "",
        "| Case | Length | FIRM | SOFT | FLOPPY | Total |",
        "|------|--------|------|------|--------|-------|",
    ]
    
    for _, row in df.iterrows():
        report_lines.append(
            f"| {row['case_name']} | {row['length']} | "
            f"{row['baseline_firm']} | {row['baseline_soft']} | "
            f"{row['baseline_floppy']} | {row['baseline_total']} |"
        )
    
    report_lines.extend([
        "",
        "## Summary",
        "",
        f"- **Mean FIRM stems:** {summary['baseline_distribution']['firm_mean']:.1f}",
        f"- **Mean SOFT stems:** {summary['baseline_distribution']['soft_mean']:.1f}",
        f"- **Mean FLOPPY stems:** {summary['baseline_distribution']['floppy_mean']:.1f}",
        "",
        "## Limitations",
        "",
        summary['note'],
        "",
        "For full robustness testing (temperature sweep, parameter sets, window jitter), ",
        "the ViennaRNA Python library would be needed to programmatically control fold ",
        "parameters. The CLI-based approach used here is limited to default parameters.",
        "",
        "**Future work:**",
        "- Temperature sweep (24°C, 37°C, 42°C)",
        "- Parameter sets (Turner 2004, Andronescu 2007, Langdon 2018)",
        "- Window boundary jitter (+/- 10-25 nt)",
        "",
    ])
    
    report_path = output_dir / "robustness_report.md"
    with open(report_path, "w") as f:
        f.write("\n".join(report_lines))
    
    print(f"✓ Saved robustness report to {report_path}")
