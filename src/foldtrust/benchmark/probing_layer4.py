"""Layer 4: SHAPE probing correlation analysis with real data."""

from pathlib import Path
from typing import Dict, List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score

from foldtrust.utils import read_fasta
from foldtrust.vienna import compute_pair_probabilities

try:
    import RNA  # noqa: F401

    HAS_RNA_PYTHON = True
except ImportError:
    HAS_RNA_PYTHON = False


def load_daslab_shape_data(shape_data_dir: Path) -> Dict[str, pd.DataFrame]:
    """
    Load SHAPE-MaP data from DasLab SARS-CoV-2 repository.

    Returns dict with 'zhang', 'incarnato', 'pyle' dataframes.
    """
    shape_files = {
        "zhang": shape_data_dir / "SHAPE data" / "zhang_invivo_reactivity.csv",
        "incarnato": shape_data_dir / "SHAPE data" / "incarnato_invivo_reactivity.csv",
        "pyle": shape_data_dir / "SHAPE data" / "pyle_reactivity.csv",
    }

    data = {}
    for name, path in shape_files.items():
        if path.exists():
            df = pd.read_csv(path)
            data[name] = df
            print(f"  Loaded {name}: {len(df)} positions")
        else:
            print(f"  Warning: {name} not found at {path}")

    return data


def extract_shape_for_region(
    shape_df: pd.DataFrame,
    start: int,
    end: int,
) -> np.ndarray:
    """
    Extract SHAPE reactivities for genomic region.

    Args:
        shape_df: DataFrame with columns including reactivity values
        start: Start position (1-indexed)
        end: End position (1-indexed, inclusive)

    Returns:
        Array of reactivity values (NaN for missing data)
    """
    region_length = end - start + 1
    reactivities = np.full(region_length, np.nan)

    # DasLab data is 1-indexed, columns are typically 'Nucleotide', 'Reactivity' or similar
    # Check actual column names
    if "Nucleotide" in shape_df.columns and "Reactivity" in shape_df.columns:
        for i, pos in enumerate(range(start, end + 1)):
            row = shape_df[shape_df["Nucleotide"] == pos]
            if len(row) > 0:
                reactivities[i] = row["Reactivity"].values[0]
    else:
        # Try to infer structure
        print(f"    Warning: Expected columns not found. Available: {list(shape_df.columns)}")
        # Assume first column is position, second is reactivity
        if len(shape_df.columns) >= 2:
            pos_col = shape_df.columns[0]
            react_col = shape_df.columns[1]
            for i, pos in enumerate(range(start, end + 1)):
                row = shape_df[shape_df[pos_col] == pos]
                if len(row) > 0:
                    reactivities[i] = row[react_col].values[0]

    return reactivities


def compute_unpaired_probabilities(prob_matrix: np.ndarray) -> np.ndarray:
    """Compute unpaired probability for each position."""
    n = prob_matrix.shape[0]
    unpaired = np.zeros(n)

    for i in range(n):
        # Sum all pairing probabilities involving position i
        paired_prob = np.sum(prob_matrix[i, :]) + np.sum(prob_matrix[:, i]) - prob_matrix[i, i]
        unpaired[i] = max(
            0.0, min(1.0, 1.0 - paired_prob / 2.0)
        )  # Divide by 2 since we counted twice

    return unpaired


def run_probing_analysis(
    output_dir: Path,
    case_dirs: List[Path],
    shape_data_dir: Path,
    figures_dir: Path,
) -> Dict:
    """
    Layer 4: SHAPE probing correlation analysis.

    Analyzes SARS-CoV-2 FSE case with DasLab SHAPE-MaP data.
    """
    print("\n### SHAPE Probing Analysis ###\n")

    if not shape_data_dir.exists():
        print(f"  Warning: SHAPE data directory not found: {shape_data_dir}")
        return {"n_cases_analyzed": 0}

    # Load SHAPE data
    shape_data = load_daslab_shape_data(shape_data_dir)

    if not shape_data:
        print("  Error: No SHAPE data loaded")
        return {"n_cases_analyzed": 0}

    results = []

    # Analyze sars2-fse case
    for case_dir in case_dirs:
        if case_dir.name != "sars2-fse":
            continue

        print(f"\nAnalyzing {case_dir.name}...")

        sequence_file = case_dir / "sequence.fa"
        _, sequence = read_fasta(sequence_file)

        # FSE is at NC_045512.2:13468-13638 (1-indexed)
        fse_start = 13468
        fse_end = 13468 + len(sequence) - 1

        print(f"  Sequence length: {len(sequence)} nt")
        print(f"  Genomic coordinates: {fse_start}-{fse_end}")

        # Compute unpaired probabilities
        prob_matrix = compute_pair_probabilities(sequence)
        unpaired_probs = compute_unpaired_probabilities(prob_matrix)

        # Extract SHAPE data for each dataset
        for dataset_name, shape_df in shape_data.items():
            print(f"\n  Dataset: {dataset_name}")

            reactivities = extract_shape_for_region(shape_df, fse_start, fse_end)

            # Count valid data points
            valid_mask = ~np.isnan(reactivities)
            n_valid = np.sum(valid_mask)

            print(f"    Valid data points: {n_valid}/{len(reactivities)}")

            if n_valid < 10:
                print("    Skipping: insufficient data")
                continue

            # Compute correlation
            rho, pval = stats.spearmanr(unpaired_probs[valid_mask], reactivities[valid_mask])

            # Compute AUROC for high reactivity
            median_react = np.nanmedian(reactivities)
            high_react = (reactivities[valid_mask] > median_react).astype(int)
            auroc = roc_auc_score(high_react, unpaired_probs[valid_mask])

            print(f"    Spearman ρ: {rho:.3f} (p={pval:.2e})")
            print(f"    AUROC: {auroc:.3f}")

            result = {
                "case": case_dir.name,
                "dataset": dataset_name,
                "length": len(sequence),
                "start": fse_start,
                "end": fse_end,
                "n_valid": n_valid,
                "spearman_rho": rho,
                "spearman_pval": pval,
                "auroc": auroc,
            }

            results.append(result)

            # Plot
            plot_shape_correlation(
                unpaired_probs,
                reactivities,
                valid_mask,
                case_dir.name,
                dataset_name,
                rho,
                pval,
                figures_dir,
            )

    # Save results
    if results:
        df = pd.DataFrame(results)
        csv_path = output_dir / "layer4_shape_probing.csv"
        df.to_csv(csv_path, index=False)
        print(f"\n✓ Saved to {csv_path}")

    summary = {
        "n_cases_analyzed": len(set([r["case"] for r in results])),
        "n_datasets": len(results),
        "datasets": list(set([r["dataset"] for r in results])),
    }

    return summary


def plot_shape_correlation(
    unpaired_probs: np.ndarray,
    reactivities: np.ndarray,
    valid_mask: np.ndarray,
    case_name: str,
    dataset_name: str,
    rho: float,
    pval: float,
    figures_dir: Path,
) -> None:
    """Plot SHAPE reactivity vs unpaired probability."""

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Scatter plot
    ax = axes[0]
    ax.scatter(
        unpaired_probs[valid_mask],
        reactivities[valid_mask],
        alpha=0.5,
        s=30,
        edgecolors="black",
        linewidth=0.5,
    )
    ax.set_xlabel("Unpaired Probability", fontsize=12, fontweight="bold")
    ax.set_ylabel("SHAPE Reactivity", fontsize=12, fontweight="bold")
    ax.set_title(
        f"{case_name} - {dataset_name}\nSpearman ρ = {rho:.3f} (p = {pval:.2e})",
        fontsize=13,
        fontweight="bold",
    )
    ax.grid(alpha=0.3, linestyle="--")

    # Position-wise plot
    ax = axes[1]
    positions = np.arange(len(unpaired_probs))

    ax.plot(
        positions,
        unpaired_probs,
        label="Unpaired probability",
        linewidth=2,
        alpha=0.8,
        color="#2E86AB",
    )

    valid_pos = positions[valid_mask]
    ax.scatter(
        valid_pos,
        reactivities[valid_mask],
        label="SHAPE reactivity",
        s=40,
        alpha=0.7,
        color="#A23B72",
        edgecolors="black",
        linewidth=0.5,
    )

    ax.set_xlabel("Position (nt)", fontsize=12, fontweight="bold")
    ax.set_ylabel("Value", fontsize=12, fontweight="bold")
    ax.set_title(f"{case_name} - {dataset_name}: Position-wise", fontsize=13, fontweight="bold")
    ax.legend(fontsize=10)
    ax.grid(alpha=0.3, linestyle="--")

    plt.tight_layout()

    filename = f"layer4_shape_{case_name}_{dataset_name}.png"
    fig_path = figures_dir / filename
    fig.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close(fig)

    print(f"    ✓ Saved figure: {filename}")
