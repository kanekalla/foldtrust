"""Robustness analysis: parameter sets, temperature, and window boundaries."""

import json
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd

from foldtrust.utils import read_fasta
from foldtrust.vienna import compute_pair_probabilities, fold_mfe, parse_stems


def run_robustness_analysis(
    output_dir: Path,
    case_paths: List[Path] = None,
    case_dirs: List[Path] = None,
    temperatures: List[int] = None,
) -> Dict:
    """
    Analyze robustness of reliability tiers under perturbations.

    Args:
        output_dir: Directory for outputs
        case_paths: List of sequence.fa file paths
        case_dirs: List of case directories with sequence.fa files (legacy)
        temperatures: List of temperatures in Celsius (default: [24, 37, 42])

    Returns:
        Dictionary with robustness results
    """
    print("\n=== Robustness Analysis ===\n")

    if temperatures is None:
        temperatures = [24, 37, 42]

    # Handle both case_paths and case_dirs for backwards compatibility
    sequence_files = []
    if case_paths:
        sequence_files = case_paths
    elif case_dirs:
        for case_dir in case_dirs:
            seq_file = case_dir / "sequence.fa"
            if seq_file.exists():
                sequence_files.append(seq_file)

    if not sequence_files:
        print("  No sequence files found")
        return {"status": "no_data"}

    results = []

    for sequence_file in sequence_files:
        case_name = sequence_file.parent.name

        print(f"Analyzing {case_name}...")

        _, sequence = read_fasta(sequence_file)

        # Baseline at 37C
        baseline_result = analyze_baseline(sequence, case_name)

        # Temperature sweep
        temp_results = analyze_temperature_sweep(sequence, case_name, temperatures)

        results.append(
            {
                "case_name": case_name,
                "length": len(sequence),
                **baseline_result,
                **temp_results,
            }
        )

    df = pd.DataFrame(results)

    csv_path = output_dir / "robustness_results.csv"
    df.to_csv(csv_path, index=False)

    summary = compute_robustness_summary(df, temperatures)

    summary_path = output_dir / "robustness_summary.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)

    generate_robustness_report(df, summary, output_dir, temperatures)

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


def analyze_temperature_sweep(sequence: str, case_name: str, temperatures: List[int]) -> Dict:
    """
    Analyze tier stability across temperatures using ViennaRNA -T flag.

    Args:
        sequence: RNA sequence
        case_name: Case identifier
        temperatures: List of temperatures in Celsius

    Returns:
        Dictionary with temperature-specific tier counts and stability metrics
    """

    baseline_temp = 37  # Standard temperature

    # Get baseline (37C) tiers
    baseline_structure, _ = fold_mfe_at_temp(sequence, baseline_temp)
    baseline_prob_matrix = compute_pair_probs_at_temp(sequence, baseline_temp)
    baseline_stems = parse_stems(baseline_structure, baseline_prob_matrix)
    baseline_tiers = {(i, j): s["flag"] for s in baseline_stems for i, j in s["pairs"]}

    results = {}

    for temp in temperatures:
        structure, mfe_energy = fold_mfe_at_temp(sequence, temp)
        prob_matrix = compute_pair_probs_at_temp(sequence, temp)
        stems = parse_stems(structure, prob_matrix)

        firm = sum(1 for s in stems if s["flag"] == "firm")
        soft = sum(1 for s in stems if s["flag"] == "soft")
        floppy = sum(1 for s in stems if s["flag"] == "floppy")

        # Compute tier stability (fraction of pairs with same tier as baseline)
        current_tiers = {(i, j): s["flag"] for s in stems for i, j in s["pairs"]}

        common_pairs = set(baseline_tiers.keys()) & set(current_tiers.keys())
        if common_pairs:
            same_tier = sum(
                1 for pair in common_pairs if baseline_tiers[pair] == current_tiers[pair]
            )
            tier_stability = same_tier / len(common_pairs)
        else:
            tier_stability = 0.0

        # Structure similarity (Jaccard index of base pairs)
        baseline_pairs = set(baseline_tiers.keys())
        current_pairs = set(current_tiers.keys())
        if baseline_pairs or current_pairs:
            structure_jaccard = len(baseline_pairs & current_pairs) / len(
                baseline_pairs | current_pairs
            )
        else:
            structure_jaccard = 1.0

        results[f"temp{temp}_firm"] = firm
        results[f"temp{temp}_soft"] = soft
        results[f"temp{temp}_floppy"] = floppy
        results[f"temp{temp}_total"] = len(stems)
        results[f"temp{temp}_tier_stability"] = tier_stability
        results[f"temp{temp}_structure_jaccard"] = structure_jaccard
        results[f"temp{temp}_mfe_energy"] = mfe_energy

    return results


def fold_mfe_at_temp(sequence: str, temperature: int) -> tuple:
    """Fold MFE structure at specified temperature."""
    import re
    import subprocess

    try:
        result = subprocess.run(
            ["RNAfold", "--noPS", "-T", str(temperature)],
            input=f">seq\n{sequence}\n",
            capture_output=True,
            text=True,
            timeout=60,
        )

        if result.returncode != 0:
            raise RuntimeError(f"RNAfold failed: {result.stderr}")

        lines = result.stdout.strip().split("\n")
        structure_line = lines[-1]

        match = re.search(r"([.()\[\]{}]+)\s+\(\s*(-?\d+\.\d+)\s*\)", structure_line)
        if not match:
            raise RuntimeError(f"Could not parse RNAfold output: {structure_line}")

        structure = match.group(1)
        energy = float(match.group(2))

        return structure, energy

    except Exception as e:
        raise RuntimeError(f"Temperature fold error: {e}")


def compute_pair_probs_at_temp(sequence: str, temperature: int) -> np.ndarray:
    """Compute pair probabilities at specified temperature."""
    import os
    import subprocess
    import tempfile

    n = len(sequence)
    prob_matrix = np.zeros((n, n))

    # Create temp directory for PS files
    with tempfile.TemporaryDirectory() as tmpdir:
        orig_dir = os.getcwd()
        os.chdir(tmpdir)

        try:
            result = subprocess.run(
                ["RNAfold", "-p", "--noPS", "-T", str(temperature)],
                input=f">seq\n{sequence}\n",
                capture_output=True,
                text=True,
                timeout=120,
            )

            if result.returncode != 0:
                raise RuntimeError(f"RNAfold partition function failed: {result.stderr}")

            dp_file = Path("seq_dp.ps")
            if not dp_file.exists():
                raise RuntimeError("RNAfold did not generate dot plot PostScript file")

            with open(dp_file, "r") as f:
                in_data_section = False
                for line in f:
                    line = line.strip()
                    if line.startswith("%start of base pair probability data"):
                        in_data_section = True
                        continue
                    elif in_data_section:
                        if line.startswith("showpage"):
                            break
                        parts = line.split()
                        if len(parts) >= 4 and parts[3] == "ubox":
                            try:
                                i = int(parts[0]) - 1
                                j = int(parts[1]) - 1
                                sqrt_prob = float(parts[2])
                                prob = sqrt_prob * sqrt_prob

                                if 0 <= i < n and 0 <= j < n:
                                    prob_matrix[i, j] = prob
                                    prob_matrix[j, i] = prob
                            except (ValueError, IndexError):
                                continue
        finally:
            os.chdir(orig_dir)

    return prob_matrix


def compute_robustness_summary(df: pd.DataFrame, temperatures: List[int]) -> Dict:
    """Compute summary statistics for robustness."""

    summary = {
        "n_cases": len(df),
        "temperatures_tested": temperatures,
        "baseline_distribution": {
            "firm_mean": float(df["baseline_firm"].mean()),
            "soft_mean": float(df["baseline_soft"].mean()),
            "floppy_mean": float(df["baseline_floppy"].mean()),
            "total_mean": float(df["baseline_total"].mean()),
        },
        "temperature_stability": {},
    }

    # Compute mean stability metrics across cases for each temperature
    for temp in temperatures:
        tier_stab_col = f"temp{temp}_tier_stability"
        struct_jacc_col = f"temp{temp}_structure_jaccard"

        if tier_stab_col in df.columns and struct_jacc_col in df.columns:
            summary["temperature_stability"][f"{temp}C"] = {
                "tier_stability_mean": float(df[tier_stab_col].mean()),
                "structure_jaccard_mean": float(df[struct_jacc_col].mean()),
            }

    return summary


def generate_robustness_report(
    df: pd.DataFrame, summary: Dict, output_dir: Path, temperatures: List[int]
) -> None:
    """Generate markdown report for robustness analysis."""

    report_lines = [
        "# Robustness Analysis",
        "",
        "## Overview",
        "",
        f"- **N cases:** {summary['n_cases']}",
        f"- **Temperatures tested:** {', '.join(str(t) + '°C' for t in temperatures)}",
        "",
        "## Baseline Tier Distribution (37°C)",
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

    report_lines.extend(
        [
            "",
            "## Temperature Stability",
            "",
            "Tier stability: fraction of base pairs maintaining the same FIRM/SOFT/FLOPPY classification.",
            "Structure Jaccard: overlap of base-paired positions between baseline and test temperature.",
            "",
            "| Case | 24°C Tier Stab. | 37°C (baseline) | 42°C Tier Stab. | 24°C Struct. Jacc. | 42°C Struct. Jacc. |",
            "|------|-----------------|-----------------|-----------------|-------------------|-------------------|",
        ]
    )

    for _, row in df.iterrows():
        temp24_tier = row.get("temp24_tier_stability", 0.0)
        temp42_tier = row.get("temp42_tier_stability", 0.0)
        temp24_jacc = row.get("temp24_structure_jaccard", 0.0)
        temp42_jacc = row.get("temp42_structure_jaccard", 0.0)

        report_lines.append(
            f"| {row['case_name']} | {temp24_tier:.3f} | 1.000 | {temp42_tier:.3f} | {temp24_jacc:.3f} | {temp42_jacc:.3f} |"
        )

    report_lines.extend(
        [
            "",
            "## Summary Statistics",
            "",
            f"- **Mean FIRM stems (37°C):** {summary['baseline_distribution']['firm_mean']:.1f}",
            f"- **Mean SOFT stems (37°C):** {summary['baseline_distribution']['soft_mean']:.1f}",
            f"- **Mean FLOPPY stems (37°C):** {summary['baseline_distribution']['floppy_mean']:.1f}",
            "",
        ]
    )

    # Temperature stability summary
    if summary.get("temperature_stability"):
        report_lines.append("### Temperature Stability (mean across cases)")
        report_lines.append("")
        report_lines.append("| Temperature | Tier Stability | Structure Jaccard |")
        report_lines.append("|-------------|----------------|-------------------|")

        for temp_label, metrics in summary["temperature_stability"].items():
            report_lines.append(
                f"| {temp_label} | {metrics['tier_stability_mean']:.3f} | {metrics['structure_jaccard_mean']:.3f} |"
            )
        report_lines.append("")

    report_lines.extend(
        [
            "",
            "## Methods",
            "",
            "- **Temperature sweep:** RNAfold with `-T` flag (24°C, 37°C, 42°C)",
            "- **Tier stability:** Fraction of common base pairs maintaining the same tier classification",
            "- **Structure Jaccard:** Jaccard index of base-paired positions",
            "",
            "## Interpretation",
            "",
            "High tier stability (close to 1.0) indicates that FIRM/SOFT/FLOPPY classifications are",
            "robust to temperature changes. High structure Jaccard indicates that the same base pairs",
            "form across temperatures. Disease windows with low stability may have competing folds",
            "sensitive to physiological temperature variation.",
            "",
        ]
    )

    report_path = output_dir / "robustness_report.md"
    with open(report_path, "w") as f:
        f.write("\n".join(report_lines))

    print(f"✓ Saved robustness report to {report_path}")
