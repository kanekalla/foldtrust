"""Calibration analysis for base-pair probabilities and reliability tiers."""

import json
from pathlib import Path
from typing import Dict, List

import matplotlib.pyplot as plt
import pandas as pd

from foldtrust.benchmark.datasets import load_reference_dataset
from foldtrust.benchmark.metrics import (
    compute_calibration_metrics,
    compute_tier_accuracy,
    parse_structure_pairs,
)
from foldtrust.vienna import compute_pair_probabilities, fold_mfe, parse_stems


def run_calibration_analysis(
    output_dir: Path,
    dataset_file: Path,
    max_sequences: int = 50,
) -> Dict:
    """
    Analyze calibration of pair probabilities and reliability tiers.

    Args:
        output_dir: Directory for outputs
        dataset_file: Path to reference dataset JSON
        max_sequences: Maximum number of sequences to analyze

    Returns:
        Dictionary with calibration results
    """
    print("\n=== Calibration Analysis ===\n")

    entries = load_reference_dataset(dataset_file)[:max_sequences]

    print(f"Analyzing calibration on {len(entries)} sequences...")

    all_bin_means = []
    all_bin_accs = []
    all_bin_counts = []

    calibration_results = []
    tier_results = []

    for i, entry in enumerate(entries):
        if (i + 1) % 10 == 0:
            print(f"  Progress: {i + 1}/{len(entries)}")

        try:
            sequence = entry["sequence"]
            ref_structure = entry["structure"]
            ref_pairs = parse_structure_pairs(ref_structure)

            prob_matrix = compute_pair_probabilities(sequence)

            calib = compute_calibration_metrics(prob_matrix, ref_pairs, n_bins=10)

            calibration_results.append(
                {
                    "name": entry["name"],
                    "length": entry["length"],
                    "ece": calib["ece"],
                    "auroc": calib["auroc"],
                    "auprc": calib["auprc"],
                }
            )

            for j, (mean, acc, count) in enumerate(
                zip(calib["bin_means"], calib["bin_accs"], calib["bin_counts"])
            ):
                all_bin_means.append(mean)
                all_bin_accs.append(acc)
                all_bin_counts.append(count)

            mfe_structure, _ = fold_mfe(sequence)
            stems = parse_stems(mfe_structure, prob_matrix)

            tier_acc = compute_tier_accuracy(stems, ref_pairs)

            tier_results.append(
                {
                    "name": entry["name"],
                    "length": entry["length"],
                    "firm_ppv": tier_acc["firm"]["ppv"],
                    "firm_count": tier_acc["firm"]["count"],
                    "soft_ppv": tier_acc["soft"]["ppv"],
                    "soft_count": tier_acc["soft"]["count"],
                    "floppy_ppv": tier_acc["floppy"]["ppv"],
                    "floppy_count": tier_acc["floppy"]["count"],
                }
            )

        except Exception as e:
            print(f"  Error processing {entry['name']}: {e}")
            continue

    print(f"\n✓ Completed {len(calibration_results)} sequences")

    calib_df = pd.DataFrame(calibration_results)
    tier_df = pd.DataFrame(tier_results)

    calib_csv = output_dir / "calibration_results.csv"
    calib_df.to_csv(calib_csv, index=False)

    tier_csv = output_dir / "tier_accuracy.csv"
    tier_df.to_csv(tier_csv, index=False)

    summary = {
        "calibration": {
            "n_sequences": len(calib_df),
            "ece_mean": float(calib_df["ece"].mean()),
            "ece_std": float(calib_df["ece"].std()),
            "auroc_mean": float(calib_df["auroc"].mean()),
            "auroc_std": float(calib_df["auroc"].std()),
            "auprc_mean": float(calib_df["auprc"].mean()),
            "auprc_std": float(calib_df["auprc"].std()),
        },
        "tier_accuracy": {
            "firm_ppv_mean": float(tier_df[tier_df["firm_count"] > 0]["firm_ppv"].mean()),
            "soft_ppv_mean": float(tier_df[tier_df["soft_count"] > 0]["soft_ppv"].mean()),
            "floppy_ppv_mean": float(tier_df[tier_df["floppy_count"] > 0]["floppy_ppv"].mean()),
            "firm_total_pairs": int(tier_df["firm_count"].sum()),
            "soft_total_pairs": int(tier_df["soft_count"].sum()),
            "floppy_total_pairs": int(tier_df["floppy_count"].sum()),
        },
    }

    summary_path = output_dir / "calibration_summary.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)

    plot_reliability_diagram(all_bin_means, all_bin_accs, all_bin_counts, output_dir)

    generate_calibration_report(summary, output_dir)

    return summary


def plot_reliability_diagram(
    bin_means: List[float], bin_accs: List[float], bin_counts: List[int], output_dir: Path
) -> None:
    """Plot reliability (calibration) diagram."""

    bin_df = pd.DataFrame(
        {
            "predicted": bin_means,
            "observed": bin_accs,
            "count": bin_counts,
        }
    )

    aggregated = (
        bin_df.groupby(pd.cut(bin_df["predicted"], bins=10))
        .agg({"predicted": "mean", "observed": "mean", "count": "sum"})
        .dropna()
    )

    fig, ax = plt.subplots(figsize=(8, 8))

    ax.plot([0, 1], [0, 1], "k--", label="Perfect calibration", linewidth=2)

    ax.scatter(
        aggregated["predicted"],
        aggregated["observed"],
        s=aggregated["count"] / 100,
        alpha=0.6,
        label="Predicted vs Observed",
    )

    ax.set_xlabel("Predicted Probability", fontsize=12)
    ax.set_ylabel("Observed Fraction Correct", fontsize=12)
    ax.set_title("Reliability Diagram (Base-Pair Probabilities)", fontsize=14, fontweight="bold")
    ax.legend()
    ax.grid(alpha=0.3)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    fig.tight_layout()

    plot_path = output_dir / "reliability_diagram.png"
    fig.savefig(plot_path, dpi=150)
    plt.close(fig)

    print(f"✓ Saved reliability diagram to {plot_path}")


def generate_calibration_report(summary: Dict, output_dir: Path) -> None:
    """Generate markdown report for calibration analysis."""

    calib = summary["calibration"]
    tier = summary["tier_accuracy"]

    report_lines = [
        "# Calibration Analysis",
        "",
        "## Base-Pair Probability Calibration",
        "",
        f"- **N sequences:** {calib['n_sequences']}",
        "",
        "| Metric | Mean | Std |",
        "|--------|------|-----|",
        f"| Expected Calibration Error (ECE) | {calib['ece_mean']:.4f} | {calib['ece_std']:.4f} |",
        f"| AUROC | {calib['auroc_mean']:.3f} | {calib['auroc_std']:.3f} |",
        f"| AUPRC | {calib['auprc_mean']:.3f} | {calib['auprc_std']:.3f} |",
        "",
        "**Interpretation:**",
        "- ECE measures the average gap between predicted probabilities and observed accuracy.",
        "- AUROC and AUPRC measure how well pair probabilities discriminate true from false pairs.",
        "",
        "## Reliability Tier Accuracy",
        "",
        "| Tier | Mean PPV | Total Pairs |",
        "|------|----------|-------------|",
        f"| **FIRM** (P ≥ 0.85) | {tier['firm_ppv_mean']:.3f} | {tier['firm_total_pairs']} |",
        f"| **SOFT** (0.5 ≤ P < 0.85) | {tier['soft_ppv_mean']:.3f} | {tier['soft_total_pairs']} |",
        f"| **FLOPPY** (P < 0.5) | {tier['floppy_ppv_mean']:.3f} | {tier['floppy_total_pairs']} |",
        "",
        "**Interpretation:**",
        "- PPV (Positive Predictive Value) shows the fraction of pairs in each tier that match the reference.",
        "- Higher PPV for FIRM tier validates FoldTrust's reliability labeling.",
        "- Lower PPV for FLOPPY tier confirms that low-probability pairs are less trustworthy.",
        "",
        "## Methods",
        "",
        "- **Calibration:** Base-pair probabilities binned into 10 intervals; observed accuracy computed per bin.",
        "- **ECE:** Sum of |predicted - observed| weighted by bin frequency.",
        "- **Tier classification:** Stems classified by mean pair probability (FIRM/SOFT/FLOPPY).",
        "",
        "See `reliability_diagram.png` for visual calibration plot.",
        "",
    ]

    report_path = output_dir / "calibration_report.md"
    with open(report_path, "w") as f:
        f.write("\n".join(report_lines))

    print(f"✓ Saved calibration report to {report_path}")
