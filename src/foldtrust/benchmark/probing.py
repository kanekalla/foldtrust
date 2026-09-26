"""Experimental probing data comparison analysis."""

import json
from pathlib import Path
from typing import Dict, Optional

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr

from foldtrust.benchmark.metrics import compute_probing_metrics
from foldtrust.vienna import compute_pair_probabilities


def run_probing_analysis(
    output_dir: Path,
    sequence: str,
    case_name: str,
    reactivity_data: Optional[Dict] = None,
) -> Dict:
    """
    Compare unpaired probabilities to experimental probing data.

    Args:
        output_dir: Directory for outputs
        sequence: RNA sequence
        case_name: Name of the case
        reactivity_data: Dict with 'reactivities' array and 'source' metadata

    Returns:
        Dictionary with correlation results
    """
    print(f"\n=== Probing Analysis: {case_name} ===\n")

    if reactivity_data is None:
        print("  No reactivity data provided; skipping probing analysis.")
        return {}

    if reactivity_data.get("source") == "MOCK_DATA_FOR_TESTING_ONLY":
        print("  Warning: Using mock data for testing only. DO NOT report these results.")
        return {}

    prob_matrix = compute_pair_probabilities(sequence)

    unpaired_probs = compute_unpaired_probabilities(prob_matrix)

    reactivities = np.array(reactivity_data["reactivities"])

    if len(reactivities) != len(sequence):
        print(
            f"  Error: Reactivity length ({len(reactivities)}) != sequence length ({len(sequence)})"
        )
        return {}

    metrics = compute_probing_metrics(unpaired_probs, reactivities)

    results = {
        "case_name": case_name,
        "sequence_length": len(sequence),
        "spearman": metrics["spearman"],
        "spearman_pvalue": metrics["spearman_pvalue"],
        "auroc": metrics["auroc"],
        "data_source": reactivity_data.get("source", "unknown"),
    }

    results_path = output_dir / f"probing_{case_name}.json"
    with open(results_path, "w") as f:
        json.dump(results, f, indent=2)

    plot_reactivity_comparison(unpaired_probs, reactivities, case_name, output_dir)

    generate_probing_report(results, output_dir, case_name)

    print(f"✓ Spearman correlation: {metrics['spearman']:.3f} (p={metrics['spearman_pvalue']:.3e})")

    return results


def compute_unpaired_probabilities(prob_matrix: np.ndarray) -> np.ndarray:
    """
    Compute unpaired probability for each nucleotide.

    Unpaired probability = 1 - sum of pairing probabilities for that position.
    """
    n = prob_matrix.shape[0]
    unpaired = np.zeros(n)

    for i in range(n):
        paired_prob = np.sum(prob_matrix[i, :]) + np.sum(prob_matrix[:, i]) - prob_matrix[i, i]
        unpaired[i] = max(0.0, 1.0 - paired_prob)

    return unpaired


def plot_reactivity_comparison(
    unpaired_probs: np.ndarray, reactivities: np.ndarray, case_name: str, output_dir: Path
) -> None:
    """Plot correlation between unpaired probabilities and experimental reactivity."""

    valid_mask = ~np.isnan(reactivities)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    ax = axes[0]
    ax.scatter(unpaired_probs[valid_mask], reactivities[valid_mask], alpha=0.5, s=20)

    spearman, p_val = spearmanr(unpaired_probs[valid_mask], reactivities[valid_mask])

    ax.set_xlabel("Unpaired Probability", fontsize=11)
    ax.set_ylabel("Experimental Reactivity (SHAPE)", fontsize=11)
    ax.set_title(
        f"{case_name}: Unpaired vs Reactivity\nSpearman ρ = {spearman:.3f} (p = {p_val:.2e})",
        fontsize=12,
        fontweight="bold",
    )
    ax.grid(alpha=0.3)

    ax = axes[1]
    positions = np.arange(len(unpaired_probs))
    ax.plot(positions, unpaired_probs, label="Unpaired probability", linewidth=1.5, alpha=0.8)
    ax.plot(
        positions[valid_mask],
        reactivities[valid_mask],
        label="Reactivity",
        linewidth=1.5,
        alpha=0.8,
    )
    ax.set_xlabel("Nucleotide Position", fontsize=11)
    ax.set_ylabel("Value", fontsize=11)
    ax.set_title(f"{case_name}: Position-wise Comparison", fontsize=12, fontweight="bold")
    ax.legend()
    ax.grid(alpha=0.3)

    fig.tight_layout()

    plot_path = output_dir / f"probing_{case_name}.png"
    fig.savefig(plot_path, dpi=150)
    plt.close(fig)

    print(f"  ✓ Saved comparison plot to {plot_path}")


def generate_probing_report(results: Dict, output_dir: Path, case_name: str) -> None:
    """Generate markdown report for probing analysis."""

    report_lines = [
        f"# Probing Data Analysis: {case_name}",
        "",
        "## Results",
        "",
        f"- **Sequence length:** {results['sequence_length']}",
        f"- **Data source:** {results['data_source']}",
        f"- **Spearman correlation:** {results['spearman']:.3f}",
        f"- **P-value:** {results['spearman_pvalue']:.2e}",
        f"- **AUROC (high reactivity):** {results['auroc']:.3f}",
        "",
        "## Interpretation",
        "",
        "- **Spearman correlation** measures monotonic relationship between unpaired probability and reactivity.",
        "- Positive correlation indicates that nucleotides predicted to be unpaired have higher reactivity.",
        "- AUROC measures discrimination of high-reactivity positions using unpaired probability.",
        "",
        "## Methods",
        "",
        "- **Unpaired probability:** 1 - (sum of pairing probabilities for each nucleotide)",
        "- **Reactivity:** Experimental SHAPE-MaP or DMS reactivity values",
        "- **AUROC:** Binary classification of high (>median) vs low reactivity",
        "",
        f"See `probing_{case_name}.png` for visualization.",
        "",
    ]

    report_path = output_dir / f"probing_{case_name}_report.md"
    with open(report_path, "w") as f:
        f.write("\n".join(report_lines))

    print(f"  ✓ Saved probing report to {report_path}")
