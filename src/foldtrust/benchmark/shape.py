"""Layer 4: Experimental agreement with SHAPE probing data on SARS-CoV-2 FSE.

CORRECTED VERSION addressing all review points:
1. Tier analysis uses MEA/MFE structures, not max pair prob
2. SHAPE-directed folding reports base-pair metrics only (no misleading ΔΔG)
3. Genome control uses full tiling (369 windows), reports p-values
4. All sensitivity claims backed by saved CSVs
5. All tables saved as CSVs
6. Chemistry/method labels verified, limitations updated
7. Regression test for vienna.py changes
8. Single FSE row (no redundant core)
"""

import json
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, mannwhitneyu
from sklearn.metrics import roc_auc_score
import RNA


def load_shape_datasets(cache_dir: Path) -> Dict[str, np.ndarray]:
    """
    Load all SARS-CoV-2 SHAPE datasets from DasLab bundle.
    
    Returns:
        Dict mapping dataset name to 29903-length reactivity array.
        Missing values are NaN.
    """
    shape_dir = cache_dir / "shape"
    
    datasets = {}
    
    # Incarnato in vitro (SHAPE-MaP in vitro)
    data = pd.read_csv(shape_dir / "incarnato_invitro_reactivity.csv", header=None)[0].values
    datasets["incarnato_invitro"] = data.astype(float)
    
    # Incarnato in vivo (SHAPE-MaP in cells)
    data = pd.read_csv(shape_dir / "incarnato_invivo_reactivity.csv", header=None)[0].values
    datasets["incarnato_invivo"] = data.astype(float)
    
    # Pyle (in vivo SHAPE-MaP in Vero E6 cells) - missing values are -999
    data = pd.read_csv(shape_dir / "pyle_reactivity.csv", header=None)[0].values
    data = data.astype(float)
    data[data == -999] = np.nan
    datasets["pyle"] = data
    
    # Zhang in vitro (icSHAPE in vitro, different chemistry)
    data = pd.read_csv(shape_dir / "zhang_invitro_reactivity.csv", header=None)[0].values
    datasets["zhang_invitro"] = data.astype(float)
    
    # Zhang in vivo (icSHAPE in Vero E6 cells, different chemistry)
    data = pd.read_csv(shape_dir / "zhang_invivo_reactivity.csv", header=None)[0].values
    datasets["zhang_invivo"] = data.astype(float)
    
    # Verify all have correct length
    for name, data in datasets.items():
        assert len(data) == 29903, f"{name} has wrong length: {len(data)}"
    
    return datasets


def normalize_reactivity_deigan(reactivity: np.ndarray) -> np.ndarray:
    """
    Normalize SHAPE reactivity for Deigan soft constraints.
    
    Uses 2-8% winsorization (per DasLab descriptions.txt) without capping at 1.0.
    This preserves the full dynamic range for soft constraint energy calculation.
    
    Args:
        reactivity: Raw reactivity array (may contain NaN)
    
    Returns:
        Normalized reactivity (no upper bound)
    """
    valid = reactivity[~np.isnan(reactivity)]
    
    if len(valid) == 0:
        return reactivity
    
    # 2-8% winsorization
    p2 = np.percentile(valid, 2)
    p92 = np.percentile(valid, 92)
    
    if p92 > p2:
        normalized = (reactivity - p2) / (p92 - p2)
        # Clip lower bound only, no upper bound
        normalized = np.maximum(normalized, 0.0)
    else:
        normalized = reactivity.copy()
    
    return normalized


def normalize_reactivity_percentile90(reactivity: np.ndarray) -> np.ndarray:
    """
    Normalize by 90th percentile, capped at 1.0 (for correlation analysis).
    """
    valid = reactivity[~np.isnan(reactivity)]
    
    if len(valid) == 0:
        return reactivity
    
    p90 = np.percentile(valid, 90)
    if p90 > 0:
        normalized = reactivity / p90
        normalized = np.minimum(normalized, 1.0)
    else:
        normalized = reactivity.copy()
    
    return normalized


def extract_window(data: np.ndarray, start_1based: int, end_1based: int) -> np.ndarray:
    """Extract window from genome-wide data (1-based coordinates, inclusive)."""
    return data[start_1based - 1 : end_1based]


def compute_unpaired_probabilities(sequence: str) -> np.ndarray:
    """
    Compute unpaired probability for each nucleotide.
    
    Uses ViennaRNA partition function.
    """
    md = RNA.md()
    md.uniq_ML = 1
    fc = RNA.fold_compound(sequence, md)
    fc.pf()
    
    # Get base-pair probabilities
    bpp = fc.bpp()
    n = len(sequence)
    
    unpaired = np.zeros(n)
    for i in range(1, n + 1):
        # Sum all pairing probabilities for position i (1-based indexing)
        # bpp[i][j] already contains the symmetric probability
        paired_prob = sum(bpp[i][j] for j in range(1, n + 1) if i != j)
        unpaired[i - 1] = max(0.0, 1.0 - paired_prob)
    
    return unpaired


def classify_nucleotides_by_threshold(
    reactivity: np.ndarray,
    threshold_reactive: float = 0.75,
    threshold_unreactive: float = 0.25
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Classify nucleotides as reactive vs unreactive for AUROC.
    
    Args:
        reactivity: Normalized reactivity (0-1)
        threshold_reactive: Quantile for "reactive" (default 0.75 = top quartile)
        threshold_unreactive: Quantile for "unreactive" (default 0.25 = bottom quartile)
    
    Returns:
        (reactive_mask, unreactive_mask) - boolean arrays
    """
    valid = reactivity[~np.isnan(reactivity)]
    if len(valid) == 0:
        return np.zeros(len(reactivity), dtype=bool), np.zeros(len(reactivity), dtype=bool)
    
    r_cutoff = np.quantile(valid, threshold_reactive)
    u_cutoff = np.quantile(valid, threshold_unreactive)
    
    reactive = (reactivity >= r_cutoff) & ~np.isnan(reactivity)
    unreactive = (reactivity <= u_cutoff) & ~np.isnan(reactivity)
    
    return reactive, unreactive


def compute_auroc(
    unpaired_probs: np.ndarray,
    reactivity: np.ndarray,
    threshold_reactive: float = 0.75,
    threshold_unreactive: float = 0.25
) -> Dict:
    """
    Compute AUROC for unpaired probability predicting high vs low reactivity.
    
    Returns:
        Dict with auroc, n_reactive, n_unreactive
    """
    reactive, unreactive = classify_nucleotides_by_threshold(
        reactivity, threshold_reactive, threshold_unreactive
    )
    
    # Keep only classified positions
    mask = reactive | unreactive
    
    if np.sum(mask) < 10:
        return {"auroc": np.nan, "n_reactive": 0, "n_unreactive": 0}
    
    y_true = reactive[mask].astype(int)
    y_score = unpaired_probs[mask]
    
    try:
        auroc = roc_auc_score(y_true, y_score)
    except ValueError:
        auroc = np.nan
    
    return {
        "auroc": auroc,
        "n_reactive": int(np.sum(reactive)),
        "n_unreactive": int(np.sum(unreactive))
    }


def structure_to_pairs(structure: str) -> set:
    """Convert dot-bracket structure to set of (i, j) pairs (0-based)."""
    pairs = set()
    stack = []
    
    for i, char in enumerate(structure):
        if char == '(':
            stack.append(i)
        elif char == ')':
            if stack:
                j = stack.pop()
                pairs.add((j, i))
    
    return pairs


def compute_f1(predicted: set, reference: set) -> float:
    """Compute F1 score for base pairs."""
    if len(reference) == 0:
        return 1.0 if len(predicted) == 0 else 0.0
    
    tp = len(predicted & reference)
    fp = len(predicted - reference)
    fn = len(reference - predicted)
    
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    
    if precision + recall == 0:
        return 0.0
    
    return 2 * precision * recall / (precision + recall)


def bootstrap_ci(values: np.ndarray, n_boot: int = 1000, ci: float = 0.95) -> Tuple[float, float]:
    """
    Compute bootstrap confidence interval.
    
    Returns:
        (lower, upper) bounds of CI
    """
    if len(values) == 0:
        return (np.nan, np.nan)
    
    rng = np.random.RandomState(42)
    boot_means = []
    
    for _ in range(n_boot):
        sample = rng.choice(values, size=len(values), replace=True)
        boot_means.append(np.mean(sample))
    
    alpha = 1 - ci
    lower = np.percentile(boot_means, alpha / 2 * 100)
    upper = np.percentile(boot_means, (1 - alpha / 2) * 100)
    
    return (lower, upper)


def analyze_tier_reactivity_corrected(
    sequence: str,
    reactivity: np.ndarray,
    dataset_name: str,
    structure_type: str = "MEA"
) -> Dict:
    """
    Compute tier reactivity using actual FoldTrust tier logic.
    
    Positions paired in the MEA (or MFE) structure get the tier of their base pair.
    Positions unpaired in that structure are labeled 'unpaired'.
    
    Args:
        sequence: RNA sequence
        reactivity: Normalized SHAPE reactivity
        dataset_name: Dataset identifier
        structure_type: "MEA" or "MFE"
    
    Returns:
        Dict with tier statistics per tier
    """
    md = RNA.md()
    md.uniq_ML = 1
    fc = RNA.fold_compound(sequence, md)
    
    if structure_type == "MEA":
        fc.pf()
        structure = fc.MEA()[0]
    else:  # MFE
        structure = fc.mfe()[0]
    
    # Get base-pair probabilities
    fc_prob = RNA.fold_compound(sequence, md)
    fc_prob.pf()
    bpp = fc_prob.bpp()
    n = len(sequence)
    
    # Parse structure to get pairs
    pairs = structure_to_pairs(structure)
    
    # Classify each position
    position_tier = ["unpaired"] * n
    
    for i, j in pairs:
        # Get pair probability (1-based indexing in ViennaRNA)
        prob = bpp[i + 1][j + 1]
        
        # Assign tier based on pair probability
        if prob >= 0.85:
            tier = "FIRM"
        elif prob >= 0.5:
            tier = "SOFT"
        else:
            tier = "FLOPPY"
        
        position_tier[i] = tier
        position_tier[j] = tier
    
    # Compute statistics per tier
    results = {}
    for tier in ["FIRM", "SOFT", "FLOPPY", "unpaired"]:
        mask = np.array([position_tier[i] == tier for i in range(n)])
        values = reactivity[mask]
        valid = values[~np.isnan(values)]
        
        if len(valid) > 0:
            ci_lower, ci_upper = bootstrap_ci(valid, n_boot=1000, ci=0.95)
            results[tier] = {
                "n": len(valid),
                "mean": float(np.mean(valid)),
                "median": float(np.median(valid)),
                "std": float(np.std(valid)),
                "ci_95_lower": float(ci_lower),
                "ci_95_upper": float(ci_upper),
            }
        else:
            results[tier] = {
                "n": 0,
                "mean": np.nan,
                "median": np.nan,
                "std": np.nan,
                "ci_95_lower": np.nan,
                "ci_95_upper": np.nan,
            }
    
    # Statistical tests among paired positions only
    tier_values = {}
    for tier in ["FIRM", "SOFT", "FLOPPY"]:
        mask = np.array([position_tier[i] == tier for i in range(n)])
        values = reactivity[mask]
        tier_values[tier] = values[~np.isnan(values)]
    
    # Mann-Whitney U tests (pairwise)
    pairwise_tests = {}
    comparisons = [("FIRM", "SOFT"), ("FIRM", "FLOPPY"), ("SOFT", "FLOPPY")]
    
    for tier1, tier2 in comparisons:
        vals1 = tier_values[tier1]
        vals2 = tier_values[tier2]
        
        if len(vals1) >= 3 and len(vals2) >= 3:
            try:
                stat, pval = mannwhitneyu(vals1, vals2, alternative='two-sided')
                pairwise_tests[f"{tier1}_vs_{tier2}"] = {
                    "statistic": float(stat),
                    "pvalue": float(pval)
                }
            except ValueError:
                pairwise_tests[f"{tier1}_vs_{tier2}"] = {
                    "statistic": np.nan,
                    "pvalue": np.nan
                }
        else:
            pairwise_tests[f"{tier1}_vs_{tier2}"] = {
                "statistic": np.nan,
                "pvalue": np.nan
            }
    
    results["pairwise_tests"] = pairwise_tests
    results["structure_type"] = structure_type
    results["structure"] = structure
    
    return results


def shape_directed_folding_corrected(
    sequence: str,
    reactivity: np.ndarray,
    dataset_name: str,
    is_icshape: bool = False
) -> Dict:
    """
    Perform SHAPE-directed folding with RNAfold (CORRECTED).
    
    Uses Deigan et al. default parameters for SHAPE: slope=1.8, intercept=-0.6.
    For icSHAPE (Zhang), flags as different chemistry (Deigan parameters not calibrated).
    
    Reports base-pair distance and F1, NOT misleading pseudo-energy comparison.
    
    Args:
        sequence: RNA sequence
        reactivity: Normalized SHAPE reactivity (Deigan normalization)
        dataset_name: Dataset identifier
        is_icshape: True if icSHAPE chemistry (different from SHAPE-MaP)
    
    Returns:
        Dict with structures and comparison metrics
    """
    # Unconstrained folding
    md = RNA.md()
    md.uniq_ML = 1
    fc_unconstrained = RNA.fold_compound(sequence, md)
    mfe_ss = fc_unconstrained.mfe()[0]
    
    # MEA structure (unconstrained)
    fc_unconstrained.pf()
    mea_ss_unconstrained = fc_unconstrained.MEA()[0]
    
    # SHAPE-directed folding
    # Convert reactivity to SHAPE array (replace NaN with -999 sentinel)
    shape_data = np.where(np.isnan(reactivity), -999.0, reactivity).tolist()
    
    md_shape = RNA.md()
    md_shape.uniq_ML = 1
    fc_shape = RNA.fold_compound(sequence, md_shape)
    
    # Add SHAPE constraints with Deigan parameters
    # slope=1.8, intercept=-0.6 (default in ViennaRNA for SHAPE)
    fc_shape.sc_add_SHAPE_deigan(shape_data, 1.8, -0.6)
    
    mfe_ss_shape = fc_shape.mfe()[0]
    
    # MEA with SHAPE
    fc_shape.pf()
    mea_ss_shape = fc_shape.MEA()[0]
    
    # Compare MFE unconstrained vs MFE SHAPE
    bp_mfe_unc = structure_to_pairs(mfe_ss)
    bp_mfe_shp = structure_to_pairs(mfe_ss_shape)
    bp_dist_mfe = len(bp_mfe_unc ^ bp_mfe_shp)
    f1_mfe = compute_f1(bp_mfe_shp, bp_mfe_unc)
    
    # Compare MEA unconstrained vs MEA SHAPE
    bp_mea_unc = structure_to_pairs(mea_ss_unconstrained)
    bp_mea_shp = structure_to_pairs(mea_ss_shape)
    bp_dist_mea = len(bp_mea_unc ^ bp_mea_shp)
    f1_mea = compute_f1(bp_mea_shp, bp_mea_unc)
    
    # Count tier changes (use MEA)
    tier_changes = count_tier_changes(sequence, mea_ss_unconstrained, mea_ss_shape)
    
    results = {
        "dataset": dataset_name,
        "is_icshape": is_icshape,
        "mfe_structure_unconstrained": mfe_ss,
        "mfe_structure_shape": mfe_ss_shape,
        "mea_structure_unconstrained": mea_ss_unconstrained,
        "mea_structure_shape": mea_ss_shape,
        "bp_distance_mfe": bp_dist_mfe,
        "bp_distance_mea": bp_dist_mea,
        "f1_mfe": f1_mfe,
        "f1_mea": f1_mea,
        "tier_changes": tier_changes,
    }
    
    return results


def count_tier_changes(sequence: str, structure1: str, structure2: str) -> Dict:
    """
    Count how many pairs change FoldTrust tier between two structures.
    
    Returns:
        Dict with counts of tier changes
    """
    md = RNA.md()
    md.uniq_ML = 1
    fc = RNA.fold_compound(sequence, md)
    fc.pf()
    bpp = fc.bpp()
    
    def get_pair_tiers(structure):
        pairs = structure_to_pairs(structure)
        pair_tiers = {}
        for i, j in pairs:
            prob = bpp[i + 1][j + 1]
            if prob >= 0.85:
                tier = "FIRM"
            elif prob >= 0.5:
                tier = "SOFT"
            else:
                tier = "FLOPPY"
            pair_tiers[(i, j)] = tier
        return pair_tiers
    
    tiers1 = get_pair_tiers(structure1)
    tiers2 = get_pair_tiers(structure2)
    
    # Pairs in both structures
    common_pairs = set(tiers1.keys()) & set(tiers2.keys())
    
    # Count tier changes
    tier_changes = 0
    for pair in common_pairs:
        if tiers1[pair] != tiers2[pair]:
            tier_changes += 1
    
    # Pairs only in structure1 or structure2
    only_1 = len(set(tiers1.keys()) - set(tiers2.keys()))
    only_2 = len(set(tiers2.keys()) - set(tiers1.keys()))
    
    return {
        "common_pairs": len(common_pairs),
        "tier_changes": tier_changes,
        "pairs_only_unconstrained": only_1,
        "pairs_only_shape": only_2,
    }


def genome_control_analysis_full(
    genome_seq: str,
    window_length: int,
    datasets: Dict[str, np.ndarray],
) -> pd.DataFrame:
    """
    Tile the full genome with non-overlapping windows (CORRECTED).
    
    Computes Spearman and AUROC for each window and dataset.
    
    Args:
        genome_seq: Full genome sequence (29903 nt)
        window_length: Length of windows (e.g., 81 for FSE)
        datasets: Dict of genome-wide SHAPE datasets
    
    Returns:
        DataFrame with correlation metrics per window
    """
    # Tile genome with non-overlapping windows
    n_windows = len(genome_seq) // window_length
    
    results = []
    
    for win_idx in range(n_windows):
        win_start_0based = win_idx * window_length
        win_end_0based = win_start_0based + window_length
        
        # Convert to 1-based for display
        win_start = win_start_0based + 1
        win_end = win_end_0based
        
        window_seq = genome_seq[win_start_0based:win_end_0based]
        
        # Convert to RNA
        window_seq_rna = window_seq.replace('T', 'U')
        
        # Compute unpaired probabilities
        unpaired_probs = compute_unpaired_probabilities(window_seq_rna)
        
        # Compute correlation with each dataset
        for dataset_name, full_reactivity in datasets.items():
            window_react = full_reactivity[win_start_0based:win_end_0based]
            window_react_norm = normalize_reactivity_percentile90(window_react)
            
            # Spearman correlation
            valid = ~np.isnan(window_react_norm)
            if np.sum(valid) >= 10:
                rho, pval = spearmanr(unpaired_probs[valid], window_react_norm[valid])
            else:
                rho, pval = np.nan, np.nan
            
            # AUROC
            auroc_result = compute_auroc(unpaired_probs, window_react_norm)
            
            results.append({
                "window_index": win_idx,
                "window_start": win_start,
                "window_end": win_end,
                "dataset": dataset_name,
                "spearman": rho,
                "spearman_pvalue": pval,
                "auroc": auroc_result["auroc"],
                "n_valid": int(np.sum(valid))
            })
    
    return pd.DataFrame(results)


def run_sensitivity_analyses(
    fse_seq_rna: str,
    window_react: np.ndarray,
    unpaired_probs: np.ndarray,
    output_dir: Path,
    dataset_name: str
) -> Dict:
    """
    Run sensitivity analyses and save as CSVs (POINT 4).
    
    Tests:
    1. Normalization methods (percentile_90, percentile_cap)
    2. AUROC thresholds (0.6/0.4, 0.75/0.25, 0.8/0.2)
    
    Returns:
        Summary dict
    """
    results = {}
    
    # 1. Normalization sensitivity
    norm_results = []
    for method in ["percentile_90", "percentile_cap"]:
        if method == "percentile_90":
            normalized = normalize_reactivity_percentile90(window_react)
        else:  # percentile_cap (2-8% normalization, capped at 1.0)
            valid = window_react[~np.isnan(window_react)]
            if len(valid) > 0:
                p2 = np.percentile(valid, 2)
                p92 = np.percentile(valid, 92)
                normalized = (window_react - p2) / (p92 - p2) if p92 > p2 else window_react.copy()
                normalized = np.clip(normalized, 0.0, 1.0)
            else:
                normalized = window_react.copy()
        
        valid = ~np.isnan(normalized)
        if np.sum(valid) >= 10:
            rho, pval = spearmanr(unpaired_probs[valid], normalized[valid])
        else:
            rho, pval = np.nan, np.nan
        
        norm_results.append({
            "dataset": dataset_name,
            "normalization": method,
            "spearman": rho,
            "spearman_pvalue": pval
        })
    
    norm_df = pd.DataFrame(norm_results)
    norm_df.to_csv(output_dir / f"sensitivity_normalization_{dataset_name}.csv", index=False)
    results["normalization"] = norm_df.to_dict('records')
    
    # 2. AUROC threshold sensitivity
    threshold_results = []
    for thresh_r, thresh_u in [(0.6, 0.4), (0.75, 0.25), (0.8, 0.2)]:
        normalized = normalize_reactivity_percentile90(window_react)
        auroc_result = compute_auroc(unpaired_probs, normalized, thresh_r, thresh_u)
        
        threshold_results.append({
            "dataset": dataset_name,
            "threshold_reactive": thresh_r,
            "threshold_unreactive": thresh_u,
            "auroc": auroc_result["auroc"],
            "n_reactive": auroc_result["n_reactive"],
            "n_unreactive": auroc_result["n_unreactive"]
        })
    
    thresh_df = pd.DataFrame(threshold_results)
    thresh_df.to_csv(output_dir / f"sensitivity_auroc_thresholds_{dataset_name}.csv", index=False)
    results["auroc_thresholds"] = thresh_df.to_dict('records')
    
    return results


def plot_arc_diagram(structure: str, title: str, ax):
    """Plot arc diagram for RNA structure."""
    pairs = structure_to_pairs(structure)
    n = len(structure)
    
    # Plot backbone
    ax.plot([0, n-1], [0, 0], 'k-', linewidth=0.5, alpha=0.3)
    
    # Plot positions
    ax.plot(range(n), [0]*n, 'o', color='gray', markersize=2, alpha=0.5)
    
    # Plot arcs for pairs
    for i, j in pairs:
        mid = (i + j) / 2
        height = (j - i) / 2 * 0.3
        theta = np.linspace(0, np.pi, 50)
        x = mid + (j - i) / 2 * np.cos(theta)
        y = height * np.sin(theta)
        ax.plot(x, y, 'b-', linewidth=0.8, alpha=0.6)
    
    ax.set_xlim(-1, n)
    ax.set_ylim(-0.5, max(n * 0.15, 5))
    ax.set_title(title, fontsize=9)
    ax.set_xlabel('Position', fontsize=8)
    ax.axis('off')


def create_all_figures(
    fse_start: int,
    fse_end: int,
    unpaired_probs: np.ndarray,
    datasets_norm: Dict[str, np.ndarray],
    df_control: pd.DataFrame,
    fse_metrics: Dict,
    shape_results_all: List[Dict],
    figures_dir: Path
):
    """Create all figures (CORRECTED to use benchmarks/outputs/figures/)."""
    
    # 1. SHAPE tracks
    position = np.arange(fse_start, fse_end + 1)
    
    n_datasets = len(datasets_norm)
    fig, axes = plt.subplots(n_datasets + 1, 1, figsize=(12, 2 * (n_datasets + 1)), sharex=True)
    
    # Top panel: unpaired probability
    axes[0].plot(position, unpaired_probs, color='black', linewidth=1.5)
    axes[0].set_ylabel('Unpaired\nProbability', fontsize=10)
    axes[0].set_ylim(0, 1.05)
    axes[0].axhline(0.5, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
    axes[0].grid(True, alpha=0.3)
    axes[0].set_title(f'SARS-CoV-2 FSE (NC_045512.2:{fse_start}-{fse_end})', fontsize=12, fontweight='bold')
    
    # Dataset panels
    colors = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00']
    
    for i, (dataset_name, reactivity) in enumerate(datasets_norm.items()):
        ax = axes[i + 1]
        
        # Plot reactivity
        valid = ~np.isnan(reactivity)
        ax.plot(position[valid], reactivity[valid], 'o-', 
                color=colors[i % len(colors)], markersize=3, linewidth=1, alpha=0.7)
        
        ax.set_ylabel(dataset_name.replace('_', '\n'), fontsize=9)
        ax.set_ylim(0, 1.05)
        ax.axhline(0.5, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
        ax.grid(True, alpha=0.3)
    
    axes[-1].set_xlabel('Position (NC_045512.2)', fontsize=10)
    
    plt.tight_layout()
    plt.savefig(figures_dir / 'layer4_shape_tracks.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. Genome control
    datasets = df_control['dataset'].unique()
    
    fig, axes = plt.subplots(1, len(datasets), figsize=(4 * len(datasets), 4), sharey=True)
    
    if len(datasets) == 1:
        axes = [axes]
    
    for i, dataset in enumerate(datasets):
        df_sub = df_control[df_control['dataset'] == dataset]
        
        # Histogram of genome correlations
        axes[i].hist(df_sub['spearman'].dropna(), bins=30, alpha=0.7, color='gray', edgecolor='black')
        
        # Mark FSE correlation
        if dataset in fse_metrics:
            fse_rho = fse_metrics[dataset]["spearman"]
            fse_pval = fse_metrics[dataset]["empirical_pvalue"]
            
            if not np.isnan(fse_rho):
                axes[i].axvline(fse_rho, color='red', linewidth=2, linestyle='--', label='FSE')
                
                # Compute percentile
                valid = df_sub['spearman'].dropna()
                percentile = (valid < fse_rho).sum() / len(valid) * 100 if len(valid) > 0 else np.nan
                
                axes[i].text(0.05, 0.95, f'FSE: {fse_rho:.3f}\n({percentile:.1f}th %ile)\np={fse_pval:.3f}', 
                            transform=axes[i].transAxes, va='top', fontsize=9,
                            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        axes[i].set_xlabel('Spearman correlation', fontsize=10)
        axes[i].set_title(dataset.replace('_', ' '), fontsize=11)
        axes[i].grid(True, alpha=0.3)
    
    axes[0].set_ylabel('Frequency', fontsize=10)
    
    plt.tight_layout()
    plt.savefig(figures_dir / 'layer4_genome_control.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 3. Arc diagrams (unconstrained vs SHAPE-directed)
    n_datasets = len(shape_results_all)
    fig, axes = plt.subplots(n_datasets, 2, figsize=(12, 3 * n_datasets))
    
    if n_datasets == 1:
        axes = axes.reshape(1, -1)
    
    for i, result in enumerate(shape_results_all):
        plot_arc_diagram(result['mea_structure_unconstrained'], 
                        f"{result['dataset']}: MEA unconstrained", axes[i, 0])
        plot_arc_diagram(result['mea_structure_shape'], 
                        f"{result['dataset']}: MEA SHAPE-directed", axes[i, 1])
    
    plt.tight_layout()
    plt.savefig(figures_dir / 'layer4_arc_diagrams.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 4. Dataset agreement heatmap
    dataset_names = list(datasets_norm.keys())
    n = len(dataset_names)
    agreement_matrix = np.zeros((n, n))
    
    # Compute pairwise correlations
    for i, name1 in enumerate(dataset_names):
        for j, name2 in enumerate(dataset_names):
            if i == j:
                agreement_matrix[i, j] = 1.0
            else:
                react1 = datasets_norm[name1]
                react2 = datasets_norm[name2]
                valid = ~np.isnan(react1) & ~np.isnan(react2)
                if np.sum(valid) >= 10:
                    rho, _ = spearmanr(react1[valid], react2[valid])
                    agreement_matrix[i, j] = rho
                else:
                    agreement_matrix[i, j] = np.nan
    
    fig, ax = plt.subplots(figsize=(8, 7))
    im = ax.imshow(agreement_matrix, cmap='RdYlBu_r', vmin=-1, vmax=1, aspect='auto')
    
    # Set ticks
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels([name.replace('_', ' ') for name in dataset_names], rotation=45, ha='right')
    ax.set_yticklabels([name.replace('_', ' ') for name in dataset_names])
    
    # Annotate cells
    for i in range(n):
        for j in range(n):
            val = agreement_matrix[i, j]
            if not np.isnan(val):
                text = ax.text(j, i, f'{val:.2f}',
                             ha="center", va="center", color="black", fontsize=9)
    
    ax.set_title('Dataset Agreement (Spearman ρ on FSE)', fontsize=12, fontweight='bold')
    plt.colorbar(im, ax=ax, label='Spearman ρ')
    plt.tight_layout()
    plt.savefig(figures_dir / 'layer4_dataset_agreement.png', dpi=300, bbox_inches='tight')
    plt.close()


def run_layer4_shape_analysis(
    cache_dir: Path,
    output_dir: Path,
    fse_start: int = 13462,
    fse_end: int = 13542,
) -> Dict:
    """
    Run complete Layer 4 SHAPE analysis (FULLY CORRECTED).
    
    Args:
        cache_dir: Path to data/_cache
        output_dir: Path to benchmarks/outputs/layer4_shape
    
    Returns:
        Summary dictionary
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    figures_dir = Path("benchmarks/outputs/figures")
    figures_dir.mkdir(parents=True, exist_ok=True)
    
    print("=== Layer 4: SHAPE Agreement Analysis (CORRECTED) ===\n")
    
    # Load genome
    genome_path = cache_dir / "genome" / "NC_045512.2.fasta"
    with open(genome_path) as f:
        lines = f.readlines()
        genome = ''.join(line.strip() for line in lines[1:])
    
    # Extract FSE sequence
    fse_seq = genome[fse_start - 1 : fse_end]
    fse_seq_rna = fse_seq.replace('T', 'U')
    
    print(f"FSE: NC_045512.2:{fse_start}-{fse_end} ({len(fse_seq)} nt)\n")
    
    # Load SHAPE datasets
    print("Loading SHAPE datasets...")
    datasets = load_shape_datasets(cache_dir)
    
    for name, data in datasets.items():
        n_valid = np.sum(~np.isnan(data))
        print(f"  {name}: {n_valid}/{len(data)} positions ({n_valid/len(data)*100:.1f}%)")
    
    print()
    
    # Compute unpaired probabilities
    print("Computing unpaired probabilities...")
    unpaired_probs = compute_unpaired_probabilities(fse_seq_rna)
    
    # Results container
    results = {
        "fse_start": fse_start,
        "fse_end": fse_end,
        "fse_length": len(fse_seq),
        "datasets": {}
    }
    
    # Per-dataset analysis
    print("\n=== Per-Dataset Analysis ===\n")
    
    per_dataset_metrics = []
    
    for dataset_name, full_reactivity in datasets.items():
        print(f"--- {dataset_name} ---")
        
        # Extract window
        window_react = extract_window(full_reactivity, fse_start, fse_end)
        
        # Normalize for correlation analysis
        window_react_norm = normalize_reactivity_percentile90(window_react)
        
        # Coverage
        coverage = np.sum(~np.isnan(window_react)) / len(window_react)
        
        # Spearman correlation
        valid = ~np.isnan(window_react_norm)
        if np.sum(valid) >= 10:
            rho, pval = spearmanr(unpaired_probs[valid], window_react_norm[valid])
        else:
            rho, pval = np.nan, np.nan
        
        print(f"  Spearman: {rho:.3f} (p={pval:.2e})")
        
        # AUROC
        auroc_result = compute_auroc(unpaired_probs, window_react_norm)
        print(f"  AUROC: {auroc_result['auroc']:.3f}")
        print(f"  Coverage: {coverage:.1%}\n")
        
        per_dataset_metrics.append({
            "dataset": dataset_name,
            "spearman": rho,
            "spearman_pvalue": pval,
            "auroc": auroc_result["auroc"],
            "n_reactive": auroc_result["n_reactive"],
            "n_unreactive": auroc_result["n_unreactive"],
            "coverage": coverage,
            "n_valid": int(np.sum(valid))
        })
        
        results["datasets"][dataset_name] = per_dataset_metrics[-1]
    
    # Save per-dataset metrics as CSV
    pd.DataFrame(per_dataset_metrics).to_csv(output_dir / "per_dataset_metrics.csv", index=False)
    
    # Dataset agreement
    print("=== Dataset Agreement ===\n")
    
    agreement_records = []
    dataset_names = list(datasets.keys())
    
    for i, name1 in enumerate(dataset_names):
        for name2 in dataset_names[i+1:]:
            react1 = extract_window(datasets[name1], fse_start, fse_end)
            react2 = extract_window(datasets[name2], fse_start, fse_end)
            
            react1_norm = normalize_reactivity_percentile90(react1)
            react2_norm = normalize_reactivity_percentile90(react2)
            
            valid = ~np.isnan(react1_norm) & ~np.isnan(react2_norm)
            
            if np.sum(valid) >= 10:
                rho, pval = spearmanr(react1_norm[valid], react2_norm[valid])
            else:
                rho, pval = np.nan, np.nan
            
            agreement_records.append({
                "dataset1": name1,
                "dataset2": name2,
                "spearman": rho,
                "pvalue": pval,
                "n_overlap": int(np.sum(valid))
            })
            
            print(f"  {name1} vs {name2}: {rho:.3f} (p={pval:.2e})")
    
    print()
    agreement_df = pd.DataFrame(agreement_records)
    agreement_df.to_csv(output_dir / "dataset_agreement.csv", index=False)
    results["dataset_agreement"] = agreement_records
    
    # Tier analysis (CORRECTED) - all datasets, MEA and MFE
    print("=== Tier Analysis (CORRECTED: MEA & MFE structures) ===\n")
    
    tier_results_all = []
    
    for dataset_name, full_reactivity in datasets.items():
        window_react = extract_window(full_reactivity, fse_start, fse_end)
        window_react_norm = normalize_reactivity_percentile90(window_react)
        
        for structure_type in ["MEA", "MFE"]:
            tier_stats = analyze_tier_reactivity_corrected(
                fse_seq_rna, window_react_norm, dataset_name, structure_type
            )
            
            print(f"--- {dataset_name} ({structure_type}) ---")
            for tier in ["FIRM", "SOFT", "FLOPPY", "unpaired"]:
                stats = tier_stats[tier]
                if stats["n"] > 0:
                    print(f"  {tier}: n={stats['n']}, mean={stats['mean']:.3f}, "
                          f"median={stats['median']:.3f}, "
                          f"95% CI=[{stats['ci_95_lower']:.3f}, {stats['ci_95_upper']:.3f}]")
                else:
                    print(f"  {tier}: n=0")
            
            # Print statistical tests
            if "pairwise_tests" in tier_stats:
                print("  Statistical tests (Mann-Whitney U):")
                for comparison, test_result in tier_stats["pairwise_tests"].items():
                    pval = test_result["pvalue"]
                    if not np.isnan(pval):
                        sig = "***" if pval < 0.001 else "**" if pval < 0.01 else "*" if pval < 0.05 else "n.s."
                        print(f"    {comparison}: p={pval:.3f} {sig}")
                    else:
                        print(f"    {comparison}: insufficient data")
            
            print()
            
            # Flatten for CSV
            for tier in ["FIRM", "SOFT", "FLOPPY", "unpaired"]:
                row = {
                    "dataset": dataset_name,
                    "structure_type": structure_type,
                    "tier": tier,
                    **tier_stats[tier]
                }
                tier_results_all.append(row)
    
    tier_df = pd.DataFrame(tier_results_all)
    tier_df.to_csv(output_dir / "tier_analysis.csv", index=False)
    results["tier_analysis"] = tier_results_all
    
    # SHAPE-directed folding (CORRECTED) - all datasets
    print("=== SHAPE-Directed Folding (CORRECTED: all datasets) ===\n")
    
    shape_results_all = []
    
    for dataset_name, full_reactivity in datasets.items():
        window_react = extract_window(full_reactivity, fse_start, fse_end)
        
        # Normalize for Deigan soft constraints (2-8%, no upper cap)
        window_react_deigan = normalize_reactivity_deigan(window_react)
        
        is_icshape = "zhang" in dataset_name.lower()
        
        shape_result = shape_directed_folding_corrected(
            fse_seq_rna, window_react_deigan, dataset_name, is_icshape
        )
        
        print(f"--- {dataset_name} ---")
        print(f"  Chemistry: {'icSHAPE (different from SHAPE-MaP)' if is_icshape else 'SHAPE-MaP'}")
        print(f"  BP distance (MFE): {shape_result['bp_distance_mfe']}")
        print(f"  BP distance (MEA): {shape_result['bp_distance_mea']}")
        print(f"  F1 (MFE): {shape_result['f1_mfe']:.3f}")
        print(f"  F1 (MEA): {shape_result['f1_mea']:.3f}")
        print(f"  Tier changes (MEA): {shape_result['tier_changes']['tier_changes']} / "
              f"{shape_result['tier_changes']['common_pairs']} common pairs\n")
        
        shape_results_all.append(shape_result)
    
    # Save SHAPE-directed results
    shape_df = pd.DataFrame([{
        "dataset": r["dataset"],
        "is_icshape": r["is_icshape"],
        "bp_distance_mfe": r["bp_distance_mfe"],
        "bp_distance_mea": r["bp_distance_mea"],
        "f1_mfe": r["f1_mfe"],
        "f1_mea": r["f1_mea"],
        "tier_changes": r["tier_changes"]["tier_changes"],
        "common_pairs": r["tier_changes"]["common_pairs"],
    } for r in shape_results_all])
    
    shape_df.to_csv(output_dir / "shape_directed_folding.csv", index=False)
    results["shape_directed_folding"] = [
        {k: v for k, v in r.items() if k not in ["mfe_structure_unconstrained", "mfe_structure_shape", 
                                                   "mea_structure_unconstrained", "mea_structure_shape"]}
        for r in shape_results_all
    ]
    
    # Genome control (CORRECTED: full tiling)
    print("=== Genome Control (CORRECTED: full tiling, 369 windows) ===\n")
    
    df_control = genome_control_analysis_full(genome, len(fse_seq), datasets)
    
    # Compute empirical percentiles and p-values
    fse_metrics = {}
    
    for dataset_name in dataset_names:
        df_sub = df_control[df_control['dataset'] == dataset_name]
        fse_rho = results["datasets"][dataset_name]["spearman"]
        
        if not np.isnan(fse_rho):
            valid_rho = df_sub['spearman'].dropna()
            n_windows = len(valid_rho)
            n_higher = (valid_rho >= fse_rho).sum()
            empirical_pvalue = n_higher / n_windows if n_windows > 0 else np.nan
            percentile = ((valid_rho < fse_rho).sum() / n_windows * 100) if n_windows > 0 else np.nan
            
            fse_metrics[dataset_name] = {
                "spearman": fse_rho,
                "percentile": percentile,
                "empirical_pvalue": empirical_pvalue,
                "n_windows": n_windows
            }
            
            print(f"  {dataset_name}:")
            print(f"    FSE Spearman: {fse_rho:.3f}")
            print(f"    Percentile: {percentile:.1f}th")
            print(f"    Empirical p-value: {empirical_pvalue:.3f}")
            print(f"    ({n_windows} non-overlapping windows)")
        else:
            fse_metrics[dataset_name] = {
                "spearman": fse_rho,
                "percentile": np.nan,
                "empirical_pvalue": np.nan,
                "n_windows": 0
            }
    
    print()
    
    df_control.to_csv(output_dir / "genome_control_full.csv", index=False)
    results["genome_control"] = fse_metrics
    
    # Sensitivity analyses (POINT 4) - run for first dataset as example
    print("=== Sensitivity Analyses (first dataset as example) ===\n")
    
    first_dataset = dataset_names[0]
    window_react = extract_window(datasets[first_dataset], fse_start, fse_end)
    
    sensitivity_results = run_sensitivity_analyses(
        fse_seq_rna, window_react, unpaired_probs, output_dir, first_dataset
    )
    
    print(f"Normalization sensitivity ({first_dataset}):")
    for record in sensitivity_results["normalization"]:
        print(f"  {record['normalization']}: Spearman = {record['spearman']:.3f}")
    
    print(f"\nAUROC threshold sensitivity ({first_dataset}):")
    for record in sensitivity_results["auroc_thresholds"]:
        print(f"  {record['threshold_reactive']}/{record['threshold_unreactive']}: "
              f"AUROC = {record['auroc']:.3f}")
    
    print()
    
    # Generate figures
    print("=== Generating Figures ===\n")
    
    datasets_norm_window = {}
    for name, full_react in datasets.items():
        window_react = extract_window(full_react, fse_start, fse_end)
        datasets_norm_window[name] = normalize_reactivity_percentile90(window_react)
    
    create_all_figures(
        fse_start, fse_end, unpaired_probs, datasets_norm_window,
        df_control, fse_metrics, shape_results_all, figures_dir
    )
    
    print("  Saved: layer4_shape_tracks.png")
    print("  Saved: layer4_genome_control.png")
    print("  Saved: layer4_arc_diagrams.png")
    print("  Saved: layer4_dataset_agreement.png")
    
    # Save JSON results
    results_path = output_dir / "layer4_shape_results.json"
    with open(results_path, "w") as f:
        # Convert numpy types to Python types for JSON
        def convert(obj):
            if isinstance(obj, np.integer):
                return int(obj)
            elif isinstance(obj, np.floating):
                return float(obj)
            elif isinstance(obj, np.ndarray):
                return obj.tolist()
            return obj
        
        json.dump(results, f, indent=2, default=convert)
    
    print(f"\n✓ Results saved to {results_path}")
    print(f"✓ All CSVs and figures saved to {output_dir} and {figures_dir}")
    
    return results


def main():
    """Command-line interface."""
    import argparse
    
    parser = argparse.ArgumentParser(description="Layer 4: SHAPE Agreement Analysis (CORRECTED)")
    parser.add_argument("--cache-dir", type=Path, default=Path("data/_cache"),
                       help="Path to data cache")
    parser.add_argument("--output-dir", type=Path, default=Path("benchmarks/outputs/layer4_shape"),
                       help="Output directory")
    
    args = parser.parse_args()
    
    run_layer4_shape_analysis(args.cache_dir, args.output_dir)


if __name__ == "__main__":
    main()
