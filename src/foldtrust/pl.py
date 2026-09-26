"""Plotting functions for FoldTrust.

Visualization tools for structure, ensemble properties, and benchmark results.
"""

from typing import Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np

from foldtrust._core import FoldData, FoldDataCollection


def arc_plot(
    fd: FoldData,
    structure_key: str = "mfe",
    color_by_tier: bool = True,
    figsize: Tuple[int, int] = (12, 4),
    save: Optional[str] = None,
    show: bool = True,
) -> plt.Figure:
    """
    Plot RNA structure as arc diagram.

    Parameters
    ----------
    fd : FoldData
        FoldData object with structure
    structure_key : str
        Which structure to plot (default: 'mfe')
    color_by_tier : bool
        Color arcs by tier (FIRM/SOFT/FLOPPY)
    figsize : tuple
        Figure size
    save : str, optional
        Path to save figure
    show : bool
        Whether to display figure

    Returns
    -------
    matplotlib.Figure
        The figure object

    Examples
    --------
    >>> fd = ft.io.read_fasta("data/cases/sars2-fse/sequence.fa")
    >>> ft.tl.run_pipeline(fd)
    >>> ft.pl.arc_plot(fd, save="sars2_arc.png")
    """
    if structure_key not in fd.structures:
        raise ValueError(f"Structure '{structure_key}' not found in fd.structures")

    structure = fd.structures[structure_key]

    # Parse base pairs
    pairs = []
    stack = []
    for i, char in enumerate(structure):
        if char == "(":
            stack.append(i)
        elif char == ")" and stack:
            j = stack.pop()
            pairs.append((j, i))

    # Get tier colors if available
    tier_colors = {"firm": "#2E7D32", "soft": "#F57C00", "floppy": "#C62828"}

    fig, ax = plt.subplots(figsize=figsize)

    # Draw sequence baseline
    positions = np.arange(len(fd.sequence))
    ax.plot(positions, np.zeros_like(positions), "k-", linewidth=0.5, alpha=0.3)

    # Draw arcs
    for i, j in pairs:
        # Determine color
        if color_by_tier and "tier" in fd.obs.columns:
            tier_i = fd.obs.loc[i, "tier"]
            tier_j = fd.obs.loc[j, "tier"]
            if tier_i == tier_j and tier_i in tier_colors:
                color = tier_colors[tier_i]
                alpha = 0.8
            else:
                color = "gray"
                alpha = 0.4
        else:
            color = "steelblue"
            alpha = 0.6

        # Draw arc
        x = np.linspace(i, j, 50)
        height = (j - i) / 2
        y = height * np.sin(np.pi * (x - i) / (j - i))
        ax.plot(x, y, color=color, alpha=alpha, linewidth=1.5)

    ax.set_xlabel("Position", fontweight="bold")
    ax.set_ylabel("Base Pairing", fontweight="bold")
    ax.set_title(f"{fd.name or 'RNA'}: {structure_key.upper()} Structure", fontweight="bold")
    ax.set_xlim(-0.5, len(fd.sequence) - 0.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Legend if coloring by tier
    if color_by_tier and "tier" in fd.obs.columns:
        from matplotlib.patches import Patch

        legend_elements = [
            Patch(facecolor=tier_colors["firm"], label="FIRM"),
            Patch(facecolor=tier_colors["soft"], label="SOFT"),
            Patch(facecolor=tier_colors["floppy"], label="FLOPPY"),
        ]
        ax.legend(handles=legend_elements, loc="upper right")

    plt.tight_layout()

    if save:
        fig.savefig(save, dpi=150, bbox_inches="tight")

    if show:
        plt.show()
    else:
        plt.close()

    return fig


def heatmap(
    fd: FoldData,
    data_key: str = "pair_probs",
    figsize: Tuple[int, int] = (8, 8),
    cmap: str = "Blues",
    save: Optional[str] = None,
    show: bool = True,
) -> plt.Figure:
    """
    Plot pair probability heatmap.

    Parameters
    ----------
    fd : FoldData
        FoldData object with pair probabilities
    data_key : str
        Which obsp data to plot (default: 'pair_probs')
    figsize : tuple
        Figure size
    cmap : str
        Colormap name
    save : str, optional
        Path to save figure
    show : bool
        Whether to display figure

    Returns
    -------
    matplotlib.Figure
        The figure object
    """
    if data_key not in fd.obsp:
        raise ValueError(f"Data '{data_key}' not found in fd.obsp")

    matrix = fd.obsp[data_key]

    fig, ax = plt.subplots(figsize=figsize)

    # Upper triangle only
    mask = np.tril(np.ones_like(matrix, dtype=bool))
    matrix_masked = np.ma.array(matrix, mask=mask)

    im = ax.imshow(matrix_masked, cmap=cmap, origin="upper", vmin=0, vmax=1)

    ax.set_xlabel("Position (j)", fontweight="bold")
    ax.set_ylabel("Position (i)", fontweight="bold")
    ax.set_title(f"{fd.name or 'RNA'}: Pair Probabilities", fontweight="bold")

    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label("Pair Probability", fontweight="bold")

    plt.tight_layout()

    if save:
        fig.savefig(save, dpi=150, bbox_inches="tight")

    if show:
        plt.show()
    else:
        plt.close()

    return fig


def tier_distribution(
    fd: FoldData,
    figsize: Tuple[int, int] = (6, 4),
    save: Optional[str] = None,
    show: bool = True,
) -> plt.Figure:
    """
    Plot tier distribution (FIRM/SOFT/FLOPPY).

    Parameters
    ----------
    fd : FoldData
        FoldData object with tier annotations
    figsize : tuple
        Figure size
    save : str, optional
        Path to save figure
    show : bool
        Whether to display figure

    Returns
    -------
    matplotlib.Figure
        The figure object
    """
    if "tier" not in fd.obs.columns:
        raise ValueError("Tier annotations not found. Run ft.tl.call_tiers first.")

    # Count tiers
    tier_counts = fd.obs["tier"].value_counts()

    tier_order = ["firm", "soft", "floppy"]
    tier_colors = {"firm": "#2E7D32", "soft": "#F57C00", "floppy": "#C62828"}

    counts = [tier_counts.get(t, 0) for t in tier_order]
    colors = [tier_colors[t] for t in tier_order]

    fig, ax = plt.subplots(figsize=figsize)

    ax.bar(tier_order, counts, color=colors, alpha=0.9, edgecolor="black", linewidth=1.5)

    ax.set_ylabel("Number of Nucleotides", fontweight="bold")
    ax.set_title(f"{fd.name or 'RNA'}: Tier Distribution", fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Add counts on bars
    for i, (t, c) in enumerate(zip(tier_order, counts)):
        if c > 0:
            ax.text(
                i, c + max(counts) * 0.02, str(int(c)), ha="center", va="bottom", fontweight="bold"
            )

    plt.tight_layout()

    if save:
        fig.savefig(save, dpi=150, bbox_inches="tight")

    if show:
        plt.show()
    else:
        plt.close()

    return fig


def shape_track(
    fd: FoldData,
    shape_key: str = "shape",
    unpaired_key: str = "unpaired_prob",
    figsize: Tuple[int, int] = (12, 5),
    save: Optional[str] = None,
    show: bool = True,
) -> plt.Figure:
    """
    Plot SHAPE reactivity and unpaired probability tracks.

    Parameters
    ----------
    fd : FoldData
        FoldData object with SHAPE data and unpaired probabilities
    shape_key : str
        Column name for SHAPE data in obs
    unpaired_key : str
        Column name for unpaired probabilities in obs
    figsize : tuple
        Figure size
    save : str, optional
        Path to save figure
    show : bool
        Whether to display figure

    Returns
    -------
    matplotlib.Figure
        The figure object
    """
    if shape_key not in fd.obs.columns:
        raise ValueError(f"SHAPE data '{shape_key}' not found in fd.obs")

    if unpaired_key not in fd.obs.columns:
        raise ValueError(f"Unpaired probabilities '{unpaired_key}' not found in fd.obs")

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize, sharex=True)

    positions = np.arange(fd.n_obs)

    # SHAPE reactivity
    shape_vals = fd.obs[shape_key].values
    valid_mask = ~np.isnan(shape_vals)

    ax1.scatter(
        positions[valid_mask],
        shape_vals[valid_mask],
        s=30,
        alpha=0.7,
        c="#A23B72",
        edgecolors="black",
        linewidth=0.5,
    )
    ax1.set_ylabel("SHAPE Reactivity", fontweight="bold")
    ax1.set_title(f"{fd.name or 'RNA'}: SHAPE vs Unpaired Probability", fontweight="bold")
    ax1.grid(alpha=0.3, linestyle="--")

    # Unpaired probability
    unpaired_vals = fd.obs[unpaired_key].values
    ax2.plot(positions, unpaired_vals, linewidth=2, color="#2E86AB", alpha=0.8)
    ax2.set_xlabel("Position", fontweight="bold")
    ax2.set_ylabel("Unpaired Probability", fontweight="bold")
    ax2.set_ylim(-0.05, 1.05)
    ax2.grid(alpha=0.3, linestyle="--")

    plt.tight_layout()

    if save:
        fig.savefig(save, dpi=150, bbox_inches="tight")

    if show:
        plt.show()
    else:
        plt.close()

    return fig


def disease_summary(
    collection: FoldDataCollection,
    figsize: Tuple[int, int] = (12, 8),
    save: Optional[str] = None,
    show: bool = True,
) -> plt.Figure:
    """
    Plot summary of tier composition across disease windows.

    Parameters
    ----------
    collection : FoldDataCollection
        Collection of FoldData objects
    figsize : tuple
        Figure size
    save : str, optional
        Path to save figure
    show : bool
        Whether to display figure

    Returns
    -------
    matplotlib.Figure
        The figure object
    """
    # Collect tier stats
    names = []
    firm_counts = []
    soft_counts = []
    floppy_counts = []

    for name, fd in collection:
        if "tier" not in fd.obs.columns:
            continue

        names.append(name)
        tier_counts = fd.obs["tier"].value_counts()
        firm_counts.append(tier_counts.get("firm", 0))
        soft_counts.append(tier_counts.get("soft", 0))
        floppy_counts.append(tier_counts.get("floppy", 0))

    if not names:
        raise ValueError("No tier annotations found in collection")

    # Plot stacked bar chart
    fig, ax = plt.subplots(figsize=figsize)

    x = np.arange(len(names))
    width = 0.6

    ax.bar(x, firm_counts, width, label="FIRM", color="#2E7D32", alpha=0.9)
    ax.bar(x, soft_counts, width, bottom=firm_counts, label="SOFT", color="#F57C00", alpha=0.9)
    ax.bar(
        x,
        floppy_counts,
        width,
        bottom=np.array(firm_counts) + np.array(soft_counts),
        label="FLOPPY",
        color="#C62828",
        alpha=0.9,
    )

    ax.set_ylabel("Number of Nucleotides", fontweight="bold", fontsize=12)
    ax.set_title("Disease Window Tier Composition", fontweight="bold", fontsize=14)
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=45, ha="right")
    ax.legend()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", alpha=0.3, linestyle="--")

    plt.tight_layout()

    if save:
        fig.savefig(save, dpi=150, bbox_inches="tight")

    if show:
        plt.show()
    else:
        plt.close()

    return fig
