"""Visualization functions for structure reports."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


def create_heatmap(prob_matrix: np.ndarray, output_path: Path):
    """
    Create heatmap visualization of base-pair probability matrix.

    Args:
        prob_matrix: NxN probability matrix
        output_path: Path to save PNG file
    """
    fig, ax = plt.subplots(figsize=(10, 9))

    mask = np.tril(np.ones_like(prob_matrix, dtype=bool), k=-1)

    display_matrix = prob_matrix.copy()
    display_matrix[mask] = np.nan

    sns.heatmap(
        display_matrix,
        cmap="YlOrRd",
        vmin=0,
        vmax=1,
        square=True,
        cbar_kws={"label": "Base-pair probability", "shrink": 0.8},
        ax=ax,
    )

    ax.set_xlabel("Position")
    ax.set_ylabel("Position")
    ax.set_title("Base-Pair Probability Matrix (upper triangle)")

    ax.invert_yaxis()

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
