"""Generate figures for Layer 5 robustness analysis."""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# Use Okabe-Ito colourblind-safe palette
OKABE_ITO = {
    "orange": "#E69F00",
    "sky_blue": "#56B4E9",
    "bluish_green": "#009E73",
    "yellow": "#F0E442",
    "blue": "#0072B2",
    "vermillion": "#D55E00",
    "reddish_purple": "#CC79A7",
    "black": "#000000",
}

plt.style.use("seaborn-v0_8-darkgrid")


def plot_temperature_retention(
    stem_retention_csv: Path,
    output_path: Path,
):
    """
    Plot stem retention heatmap across temperatures.

    Args:
        stem_retention_csv: Path to temperature_stem_retention.csv
        output_path: Output PNG path
    """
    df = pd.read_csv(stem_retention_csv)

    # Pivot for heatmap: rows=case+tier, cols=temperature
    df["case_tier"] = df["case"] + " " + df["tier"]
    pivot = df.pivot_table(
        index="case_tier", columns="temperature", values="retention", aggfunc="mean"
    )

    # Sort by case order
    case_order = ["sars2-fse", "smn2-iss-n1", "cftr-5utr", "mapt-e10", "hcv-ires-dii"]
    tier_order = ["FIRM", "SOFT", "FLOPPY"]

    ordered_rows = []
    for case in case_order:
        for tier in tier_order:
            row_name = f"{case} {tier}"
            if row_name in pivot.index:
                ordered_rows.append(row_name)

    pivot = pivot.reindex(ordered_rows)

    # Plot heatmap
    fig, ax = plt.subplots(figsize=(8, 10))

    sns.heatmap(
        pivot,
        annot=True,
        fmt=".2f",
        cmap="RdYlGn",
        vmin=0.0,
        vmax=1.0,
        cbar_kws={"label": "Retention"},
        ax=ax,
    )

    ax.set_xlabel("Temperature", fontsize=12)
    ax.set_ylabel("Case - Tier", fontsize=12)
    ax.set_title(
        "Stem Retention Across Temperatures\n(vs 37°C Turner2004 baseline)",
        fontsize=14,
        fontweight="bold",
    )

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"✓ {output_path.name}")


def plot_params_retention(
    stem_retention_csv: Path,
    output_path: Path,
):
    """
    Plot stem retention heatmap across parameter sets.

    Args:
        stem_retention_csv: Path to parameters_stem_retention.csv
        output_path: Output PNG path
    """
    df = pd.read_csv(stem_retention_csv)

    # Pivot for heatmap
    df["case_tier"] = df["case"] + " " + df["tier"]
    pivot = df.pivot_table(
        index="case_tier", columns="parameters", values="retention", aggfunc="mean"
    )

    # Sort by case order
    case_order = ["sars2-fse", "smn2-iss-n1", "cftr-5utr", "mapt-e10", "hcv-ires-dii"]
    tier_order = ["FIRM", "SOFT", "FLOPPY"]

    ordered_rows = []
    for case in case_order:
        for tier in tier_order:
            row_name = f"{case} {tier}"
            if row_name in pivot.index:
                ordered_rows.append(row_name)

    pivot = pivot.reindex(ordered_rows)

    # Plot heatmap
    fig, ax = plt.subplots(figsize=(8, 10))

    sns.heatmap(
        pivot,
        annot=True,
        fmt=".2f",
        cmap="RdYlGn",
        vmin=0.0,
        vmax=1.0,
        cbar_kws={"label": "Retention"},
        ax=ax,
    )

    ax.set_xlabel("Parameter Set", fontsize=12)
    ax.set_ylabel("Case - Tier", fontsize=12)
    ax.set_title(
        "Stem Retention Across Parameter Sets\n(vs Turner2004 baseline, 37°C)",
        fontsize=14,
        fontweight="bold",
    )

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"✓ {output_path.name}")


def plot_window_retention(
    stem_retention_csv: Path,
    output_path: Path,
):
    """
    Plot core-stem retention vs flank size.

    Args:
        stem_retention_csv: Path to window_context_stem_retention.csv
        output_path: Output PNG path
    """
    df = pd.read_csv(stem_retention_csv)

    # Calculate mean retention per tier per flank size
    summary = df.groupby(["flank_size", "tier"])["retention"].mean().reset_index()

    # Plot
    fig, ax = plt.subplots(figsize=(10, 6))

    colors = {
        "FIRM": OKABE_ITO["blue"],
        "SOFT": OKABE_ITO["orange"],
        "FLOPPY": OKABE_ITO["bluish_green"],
    }

    for tier in ["FIRM", "SOFT", "FLOPPY"]:
        tier_data = summary[summary["tier"] == tier]
        if len(tier_data) > 0:
            ax.plot(
                tier_data["flank_size"],
                tier_data["retention"],
                marker="o",
                linewidth=2.5,
                markersize=8,
                label=tier,
                color=colors.get(tier, "gray"),
            )

    ax.set_xlabel("Flank Size (nt on each side)", fontsize=12)
    ax.set_ylabel("Core-Window Stem Retention", fontsize=12)
    ax.set_title(
        "Effect of Window Context on Stem Retention\n(Turner2004, 37°C)",
        fontsize=14,
        fontweight="bold",
    )
    ax.set_ylim(0.0, 1.05)
    ax.legend(title="Tier", fontsize=10, title_fontsize=11)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"✓ {output_path.name}")


def generate_all_layer5_figures(
    output_dir: Path,
    figures_dir: Path,
):
    """
    Generate all Layer 5 figures.

    Args:
        output_dir: benchmarks/outputs/layer5
        figures_dir: benchmarks/outputs/figures
    """
    figures_dir.mkdir(parents=True, exist_ok=True)

    print("\nGenerating Layer 5 figures...")

    # Temperature retention
    temp_csv = output_dir / "temperature_stem_retention.csv"
    if temp_csv.exists():
        plot_temperature_retention(temp_csv, figures_dir / "layer5_temperature_retention.png")

    # Parameter set retention
    params_csv = output_dir / "parameters_stem_retention.csv"
    if params_csv.exists():
        plot_params_retention(params_csv, figures_dir / "layer5_params_retention.png")

    # Window context retention
    window_csv = output_dir / "window_context_stem_retention.csv"
    if window_csv.exists():
        plot_window_retention(window_csv, figures_dir / "layer5_window_retention.png")

    print("✓ All Layer 5 figures generated")


if __name__ == "__main__":
    generate_all_layer5_figures(
        Path("benchmarks/outputs/layer5"), Path("benchmarks/outputs/figures")
    )
