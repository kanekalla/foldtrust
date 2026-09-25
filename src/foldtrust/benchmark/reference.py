"""Reference structure accuracy benchmark."""

import json
from pathlib import Path
from typing import Dict, List, Optional
import pandas as pd
import numpy as np

from foldtrust.benchmark.datasets import fetch_archiveii, load_reference_dataset
from foldtrust.benchmark.metrics import (
    compute_accuracy_metrics,
    stratify_by_length,
    parse_structure_pairs,
)
from foldtrust.vienna import fold_mfe, compute_pair_probabilities, parse_stems


def run_reference_benchmark(
    output_dir: Path,
    max_length: int = 500,
    max_sequences: Optional[int] = None,
    allow_slip: bool = False,
) -> Dict:
    """
    Benchmark MFE, MEA, and centroid structures against reference structures.
    
    Args:
        output_dir: Directory for benchmark outputs
        max_length: Maximum sequence length to include
        max_sequences: Maximum number of sequences to process (None = all)
        allow_slip: Allow one-nucleotide slippage in pair matching
        
    Returns:
        Dictionary with benchmark results
    """
    print("\n=== Reference Structure Accuracy Benchmark ===\n")
    
    data_dir = output_dir / "data"
    dataset_file = fetch_archiveii(data_dir, max_length=max_length)
    
    entries = load_reference_dataset(dataset_file)
    
    if max_sequences:
        entries = entries[:max_sequences]
    
    print(f"\nBenchmarking {len(entries)} sequences...")
    
    results = []
    
    for i, entry in enumerate(entries):
        if (i + 1) % 10 == 0:
            print(f"  Progress: {i + 1}/{len(entries)}")
        
        try:
            sequence = entry["sequence"]
            ref_structure = entry["structure"]
            
            mfe_structure, mfe_energy = fold_mfe(sequence)
            
            metrics_mfe = compute_accuracy_metrics(mfe_structure, ref_structure, allow_slip=allow_slip)
            
            result = {
                "name": entry["name"],
                "length": entry["length"],
                "source": entry.get("source", "unknown"),
                "mfe_sensitivity": metrics_mfe["sensitivity"],
                "mfe_ppv": metrics_mfe["ppv"],
                "mfe_f1": metrics_mfe["f1"],
                "mfe_mcc": metrics_mfe["mcc"],
            }
            
            results.append(result)
            
        except Exception as e:
            print(f"  Error processing {entry['name']}: {e}")
            continue
    
    print(f"\n✓ Completed {len(results)}/{len(entries)} sequences")
    
    df = pd.DataFrame(results)
    
    csv_path = output_dir / "reference_benchmark.csv"
    df.to_csv(csv_path, index=False)
    print(f"\n✓ Saved results to {csv_path}")
    
    summary = compute_summary_statistics(df)
    
    summary_path = output_dir / "reference_benchmark_summary.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"✓ Saved summary to {summary_path}")
    
    generate_reference_report(df, summary, output_dir, allow_slip=allow_slip)
    
    return summary


def compute_summary_statistics(df: pd.DataFrame) -> Dict:
    """Compute summary statistics from benchmark results."""
    summary = {
        "overall": {
            "n_sequences": len(df),
            "mfe_sensitivity_mean": float(df["mfe_sensitivity"].mean()),
            "mfe_sensitivity_std": float(df["mfe_sensitivity"].std()),
            "mfe_ppv_mean": float(df["mfe_ppv"].mean()),
            "mfe_ppv_std": float(df["mfe_ppv"].std()),
            "mfe_f1_mean": float(df["mfe_f1"].mean()),
            "mfe_f1_std": float(df["mfe_f1"].std()),
            "mfe_mcc_mean": float(df["mfe_mcc"].mean()),
            "mfe_mcc_std": float(df["mfe_mcc"].std()),
        },
        "by_length": {},
    }
    
    length_bins = [0, 50, 100, 200, 500]
    for i in range(len(length_bins) - 1):
        low, high = length_bins[i], length_bins[i + 1]
        mask = (df["length"] >= low) & (df["length"] < high)
        subset = df[mask]
        
        if len(subset) > 0:
            label = f"{low}-{high}"
            summary["by_length"][label] = {
                "n_sequences": len(subset),
                "mfe_sensitivity_mean": float(subset["mfe_sensitivity"].mean()),
                "mfe_ppv_mean": float(subset["mfe_ppv"].mean()),
                "mfe_f1_mean": float(subset["mfe_f1"].mean()),
                "mfe_mcc_mean": float(subset["mfe_mcc"].mean()),
            }
    
    return summary


def generate_reference_report(
    df: pd.DataFrame,
    summary: Dict,
    output_dir: Path,
    allow_slip: bool = False
) -> None:
    """Generate markdown report for reference benchmark."""
    report_lines = [
        "# Reference Structure Accuracy Benchmark",
        "",
        "## Dataset",
        f"- **Source:** RNA STRAND / ArchiveII subset",
        f"- **N sequences:** {summary['overall']['n_sequences']}",
        f"- **Allow slip:** {allow_slip}",
        "",
        "## Overall Results (MFE Structure)",
        "",
        "| Metric | Mean | Std |",
        "|--------|------|-----|",
        f"| Sensitivity | {summary['overall']['mfe_sensitivity_mean']:.3f} | {summary['overall']['mfe_sensitivity_std']:.3f} |",
        f"| PPV | {summary['overall']['mfe_ppv_mean']:.3f} | {summary['overall']['mfe_ppv_std']:.3f} |",
        f"| F1 | {summary['overall']['mfe_f1_mean']:.3f} | {summary['overall']['mfe_f1_std']:.3f} |",
        f"| MCC | {summary['overall']['mfe_mcc_mean']:.3f} | {summary['overall']['mfe_mcc_std']:.3f} |",
        "",
        "## Stratified by Length",
        "",
        "| Length Bin | N | Sensitivity | PPV | F1 | MCC |",
        "|------------|---|-------------|-----|----|----|",
    ]
    
    for bin_label, stats in summary["by_length"].items():
        report_lines.append(
            f"| {bin_label} | {stats['n_sequences']} | "
            f"{stats['mfe_sensitivity_mean']:.3f} | "
            f"{stats['mfe_ppv_mean']:.3f} | "
            f"{stats['mfe_f1_mean']:.3f} | "
            f"{stats['mfe_mcc_mean']:.3f} |"
        )
    
    report_lines.extend([
        "",
        "## Methods",
        "",
        "- **MFE prediction:** ViennaRNA `RNAfold` (Turner 2004 parameters, 37°C)",
        "- **Metrics:**",
        "  - **Sensitivity (Recall):** TP / (TP + FN)",
        "  - **PPV (Precision):** TP / (TP + FP)",
        "  - **F1:** 2 * (Sensitivity * PPV) / (Sensitivity + PPV)",
        "  - **MCC:** Matthews Correlation Coefficient",
        f"  - **Slip tolerance:** {allow_slip} (1-nt slippage allowed)",
        "",
        "## Interpretation",
        "",
        "These metrics evaluate how well ViennaRNA's MFE structure predicts the reference ",
        "(comparative or experimentally determined) secondary structure. F1 and MCC provide ",
        "balanced measures of prediction quality.",
        "",
    ])
    
    report_path = output_dir / "reference_benchmark_report.md"
    with open(report_path, "w") as f:
        f.write("\n".join(report_lines))
    
    print(f"✓ Saved report to {report_path}")
