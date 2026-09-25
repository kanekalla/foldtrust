"""Layer 4: Experimental agreement with SHAPE probing data on SARS-CoV-2 FSE."""

import json
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score, roc_curve
import RNA

from foldtrust.vienna import compute_pair_probabilities


def load_shape_datasets(cache_dir: Path) -> Dict[str, np.ndarray]:
    """
    Load all SARS-CoV-2 SHAPE datasets from DasLab bundle.
    
    Returns:
        Dict mapping dataset name to 29903-length reactivity array.
        Missing values are NaN.
    """
    shape_dir = cache_dir / "shape"
    
    datasets = {}
    
    # Incarnato in vitro
    data = pd.read_csv(shape_dir / "incarnato_invitro_reactivity.csv", header=None)[0].values
    datasets["incarnato_invitro"] = data.astype(float)
    
    # Incarnato in vivo
    data = pd.read_csv(shape_dir / "incarnato_invivo_reactivity.csv", header=None)[0].values
    datasets["incarnato_invivo"] = data.astype(float)
    
    # Pyle (in vivo SHAPE-MaP) - missing values are -999
    data = pd.read_csv(shape_dir / "pyle_reactivity.csv", header=None)[0].values
    data = data.astype(float)
    data[data == -999] = np.nan
    datasets["pyle"] = data
    
    # Zhang in vitro (icSHAPE)
    data = pd.read_csv(shape_dir / "zhang_invitro_reactivity.csv", header=None)[0].values
    datasets["zhang_invitro"] = data.astype(float)
    
    # Zhang in vivo (icSHAPE)
    data = pd.read_csv(shape_dir / "zhang_invivo_reactivity.csv", header=None)[0].values
    datasets["zhang_invivo"] = data.astype(float)
    
    # Verify all have correct length
    for name, data in datasets.items():
        assert len(data) == 29903, f"{name} has wrong length: {len(data)}"
    
    return datasets


def normalize_reactivity(reactivity: np.ndarray, method: str = "percentile_90") -> np.ndarray:
    """
    Normalize SHAPE reactivity data.
    
    Args:
        reactivity: Raw reactivity array (may contain NaN)
        method: Normalization method
            - "percentile_90": Scale by 90th percentile (cap at 1.0)
            - "percentile_cap": 2-8% winsorization (Incarnato scale)
    
    Returns:
        Normalized reactivity (0-1 scale)
    """
    valid = reactivity[~np.isnan(reactivity)]
    
    if len(valid) == 0:
        return reactivity
    
    if method == "percentile_90":
        p90 = np.percentile(valid, 90)
        if p90 > 0:
            normalized = reactivity / p90
            normalized = np.minimum(normalized, 1.0)
        else:
            normalized = reactivity.copy()
    
    elif method == "percentile_cap":
        # 2-8% winsorization as described in DasLab descriptions.txt
        p2 = np.percentile(valid, 2)
        p92 = np.percentile(valid, 92)
        
        normalized = (reactivity - p2) / (p92 - p2) if p92 > p2 else reactivity.copy()
        normalized = np.clip(normalized, 0.0, 1.0)
    
    else:
        raise ValueError(f"Unknown normalization method: {method}")
    
    return normalized


def extract_window(data: np.ndarray, start_1based: int, end_1based: int) -> np.ndarray:
    """Extract window from genome-wide data (1-based coordinates, inclusive)."""
    return data[start_1based - 1 : end_1based]


def compute_unpaired_probabilities(sequence: str) -> np.ndarray:
    """
    Compute unpaired probability for each nucleotide.
    
    Uses ViennaRNA partition function.
    """
    prob_matrix = compute_pair_probabilities(sequence)
    n = len(sequence)
    unpaired = np.zeros(n)
    
    for i in range(n):
        # Sum all pairing probabilities involving position i
        paired_prob = np.sum(prob_matrix[i, :]) + np.sum(prob_matrix[:, i]) - prob_matrix[i, i]
        unpaired[i] = max(0.0, 1.0 - paired_prob)
    
    return unpaired


def classify_nucleotides_by_threshold(
    reactivity: np.ndarray,
    threshold_reactive: float = 0.75,
    threshold_unreactive: float = 0.25
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Classify nucleotides as reactive vs unreactive.
    
    Args:
        reactivity: Normalized reactivity (0-1)
        threshold_reactive: Quantile for "reactive" (e.g., 0.75 = top quartile)
        threshold_unreactive: Quantile for "unreactive" (e.g., 0.25 = bottom quartile)
    
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
    
    if np.sum(mask) < 10:  # Need enough positions
        return {"auroc": np.nan, "n_reactive": 0, "n_unreactive": 0}
    
    y_true = reactive[mask].astype(int)
    y_score = unpaired_probs[mask]
    
    try:
        auroc = roc_auc_score(y_true, y_score)
    except ValueError:
        auroc = np.nan
    
    return {
        "auroc": auroc,
        "n_reactive": np.sum(reactive),
        "n_unreactive": np.sum(unreactive)
    }


def analyze_tier_reactivity(
    unpaired_probs: np.ndarray,
    pair_probs: np.ndarray,
    reactivity: np.ndarray,
    tier_thresholds: Dict[str, Tuple[float, float]] = None
) -> Dict:
    """
    Compute mean reactivity for nucleotides in FIRM/SOFT/FLOPPY pairs vs unpaired.
    
    Args:
        unpaired_probs: Unpaired probability per position
        pair_probs: Pair probability matrix
        reactivity: SHAPE reactivity per position
        tier_thresholds: Dict with tier names -> (min_prob, max_prob)
    
    Returns:
        Dict with mean reactivity per tier
    """
    if tier_thresholds is None:
        tier_thresholds = {
            "FIRM": (0.85, 1.0),
            "SOFT": (0.5, 0.85),
            "FLOPPY": (0.0, 0.5)
        }
    
    n = len(unpaired_probs)
    
    # Classify each position by its strongest pairing tier
    position_tier = ["unpaired"] * n
    max_pair_prob = np.zeros(n)
    
    for i in range(n):
        for j in range(n):
            if i != j and pair_probs[i, j] > max_pair_prob[i]:
                max_pair_prob[i] = pair_probs[i, j]
    
    for i in range(n):
        if max_pair_prob[i] >= 0.85:
            position_tier[i] = "FIRM"
        elif max_pair_prob[i] >= 0.5:
            position_tier[i] = "SOFT"
        elif max_pair_prob[i] > 0:
            position_tier[i] = "FLOPPY"
        # else stays "unpaired"
    
    # Compute mean reactivity per tier
    results = {}
    for tier in ["FIRM", "SOFT", "FLOPPY", "unpaired"]:
        mask = [position_tier[i] == tier for i in range(n)]
        values = reactivity[mask]
        valid = values[~np.isnan(values)]
        
        results[tier] = {
            "mean_reactivity": np.mean(valid) if len(valid) > 0 else np.nan,
            "std_reactivity": np.std(valid) if len(valid) > 0 else np.nan,
            "n_positions": len(valid)
        }
    
    return results


def shape_directed_folding(
    sequence: str,
    reactivity: np.ndarray,
    reference_structure: str = None
) -> Dict:
    """
    Perform SHAPE-directed folding with RNAfold.
    
    Uses Deigan et al. default parameters: slope=1.8, intercept=-0.6
    
    Args:
        sequence: RNA sequence
        reactivity: Normalized SHAPE reactivity (0-1)
        reference_structure: Optional reference structure for comparison
    
    Returns:
        Dict with structures and comparison metrics
    """
    # Unconstrained folding
    md = RNA.md()
    md.uniq_ML = 1  # Unique multiloop decomposition
    fc_unconstrained = RNA.fold_compound(sequence, md)
    mfe_result = fc_unconstrained.mfe()
    mfe_ss = mfe_result[0]
    mfe_unconstrained = mfe_result[1] if isinstance(mfe_result[1], float) else mfe_result[0]
    mfe_unconstrained = mfe_result[1] if len(mfe_result) > 1 else 0.0
    
    # MEA structure (unconstrained)
    fc_unconstrained.pf()
    mea_result = fc_unconstrained.MEA()
    mea_ss_unconstrained = mea_result[0]
    mea_energy = mea_result[1] if len(mea_result) > 1 else 0.0
    
    # SHAPE-directed folding
    # Convert reactivity to SHAPE array (replace NaN with -999 sentinel)
    shape_data = np.where(np.isnan(reactivity), -999.0, reactivity).tolist()
    
    md_shape = RNA.md()
    md_shape.uniq_ML = 1
    fc_shape = RNA.fold_compound(sequence, md_shape)
    
    # Add SHAPE constraints with Deigan parameters
    # slope=1.8, intercept=-0.6 (default in ViennaRNA for SHAPE)
    fc_shape.sc_add_SHAPE_deigan(shape_data, 1.8, -0.6)
    
    mfe_shape_result = fc_shape.mfe()
    mfe_ss_shape = mfe_shape_result[0]
    mfe_shape = mfe_shape_result[1] if len(mfe_shape_result) > 1 else 0.0
    
    # MEA with SHAPE
    fc_shape.pf()
    mea_shape_result = fc_shape.MEA()
    mea_ss_shape = mea_shape_result[0]
    mea_energy_shape = mea_shape_result[1] if len(mea_shape_result) > 1 else 0.0
    
    # Compute base-pair distance and F1 if reference provided
    results = {
        "mfe_unconstrained": mfe_unconstrained,
        "mfe_structure": mfe_ss,
        "mea_structure_unconstrained": mea_ss_unconstrained,
        "mfe_shape": mfe_shape,
        "mfe_structure_shape": mfe_ss_shape,
        "mea_structure_shape": mea_ss_shape,
    }
    
    # Compare MFE unconstrained vs MFE SHAPE (compute manually)
    bp_mfe_unc = structure_to_pairs(mfe_ss)
    bp_mfe_shp = structure_to_pairs(mfe_ss_shape)
    bp_dist_mfe = len(bp_mfe_unc ^ bp_mfe_shp)  # Symmetric difference
    results["bp_distance_mfe"] = bp_dist_mfe
    
    # Compare MEA unconstrained vs MEA SHAPE
    bp_mea_unc = structure_to_pairs(mea_ss_unconstrained)
    bp_mea_shp = structure_to_pairs(mea_ss_shape)
    bp_dist_mea = len(bp_mea_unc ^ bp_mea_shp)
    results["bp_distance_mea"] = bp_dist_mea
    
    # Compute F1 if reference provided
    if reference_structure:
        # MFE vs reference
        bp_mfe_unc = structure_to_pairs(mfe_ss)
        bp_mfe_shape = structure_to_pairs(mfe_ss_shape)
        bp_ref = structure_to_pairs(reference_structure)
        
        results["f1_mfe_unconstrained_vs_ref"] = compute_f1(bp_mfe_unc, bp_ref)
        results["f1_mfe_shape_vs_ref"] = compute_f1(bp_mfe_shape, bp_ref)
    
    return results


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


def genome_control_analysis(
    genome_seq: str,
    window_length: int,
    n_windows: int,
    datasets: Dict[str, np.ndarray],
    seed: int = 42
) -> pd.DataFrame:
    """
    Sample random windows across genome and compute SHAPE-structure correlation.
    
    Args:
        genome_seq: Full genome sequence (29903 nt)
        window_length: Length of windows (e.g., 81 for FSE)
        n_windows: Number of random windows to sample
        datasets: Dict of genome-wide SHAPE datasets
        seed: Random seed
    
    Returns:
        DataFrame with correlation metrics per window
    """
    np.random.seed(seed)
    
    # Valid window positions (1-based start, ensuring full window fits)
    max_start = len(genome_seq) - window_length + 1
    window_starts = np.random.choice(max_start, size=n_windows, replace=False) + 1  # Convert to 1-based
    
    results = []
    
    for win_start in window_starts:
        win_end = win_start + window_length - 1
        window_seq = genome_seq[win_start - 1 : win_end]
        
        # Convert to RNA
        window_seq_rna = window_seq.replace('T', 'U')
        
        # Compute unpaired probabilities
        unpaired_probs = compute_unpaired_probabilities(window_seq_rna)
        
        # Compute correlation with each dataset
        for dataset_name, full_reactivity in datasets.items():
            window_react = extract_window(full_reactivity, win_start, win_end)
            window_react_norm = normalize_reactivity(window_react, method="percentile_90")
            
            # Spearman correlation
            valid = ~np.isnan(window_react_norm)
            if np.sum(valid) >= 10:
                rho, pval = spearmanr(unpaired_probs[valid], window_react_norm[valid])
            else:
                rho, pval = np.nan, np.nan
            
            results.append({
                "window_start": win_start,
                "window_end": win_end,
                "dataset": dataset_name,
                "spearman": rho,
                "spearman_pvalue": pval,
                "n_valid": np.sum(valid)
            })
    
    return pd.DataFrame(results)


def plot_shape_tracks(
    position: np.ndarray,
    unpaired_probs: np.ndarray,
    datasets_norm: Dict[str, np.ndarray],
    output_path: Path,
    title: str = "SARS-CoV-2 FSE"
):
    """
    Plot SHAPE reactivity tracks vs unpaired probability.
    
    One panel per dataset.
    """
    n_datasets = len(datasets_norm)
    
    fig, axes = plt.subplots(n_datasets + 1, 1, figsize=(12, 2 * (n_datasets + 1)), sharex=True)
    
    # Top panel: unpaired probability
    axes[0].plot(position, unpaired_probs, color='black', linewidth=1.5)
    axes[0].set_ylabel('Unpaired\nProbability', fontsize=10)
    axes[0].set_ylim(0, 1.05)
    axes[0].axhline(0.5, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
    axes[0].grid(True, alpha=0.3)
    axes[0].set_title(title, fontsize=12, fontweight='bold')
    
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
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def plot_genome_control(
    df_control: pd.DataFrame,
    fse_spearman: Dict[str, float],
    output_path: Path
):
    """
    Plot genome-control distribution with FSE marked.
    """
    datasets = df_control['dataset'].unique()
    
    fig, axes = plt.subplots(1, len(datasets), figsize=(4 * len(datasets), 4), sharey=True)
    
    if len(datasets) == 1:
        axes = [axes]
    
    for i, dataset in enumerate(datasets):
        df_sub = df_control[df_control['dataset'] == dataset]
        
        # Histogram of genome correlations
        axes[i].hist(df_sub['spearman'].dropna(), bins=20, alpha=0.7, color='gray', edgecolor='black')
        
        # Mark FSE correlation
        if dataset in fse_spearman:
            fse_rho = fse_spearman[dataset]
            axes[i].axvline(fse_rho, color='red', linewidth=2, linestyle='--', label='FSE')
            
            # Compute percentile
            valid = df_sub['spearman'].dropna()
            percentile = (valid < fse_rho).sum() / len(valid) * 100 if len(valid) > 0 else np.nan
            axes[i].text(0.05, 0.95, f'FSE: {fse_rho:.3f}\n({percentile:.1f}th %ile)', 
                        transform=axes[i].transAxes, va='top', fontsize=9,
                        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        axes[i].set_xlabel('Spearman correlation', fontsize=10)
        axes[i].set_title(dataset.replace('_', ' '), fontsize=11)
        axes[i].grid(True, alpha=0.3)
    
    axes[0].set_ylabel('Frequency', fontsize=10)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def run_layer4_shape_analysis(
    cache_dir: Path,
    output_dir: Path,
    fse_start: int = 13462,
    fse_end: int = 13542,
    core_start: int = 13462,
    core_end: int = 13542,
    n_genome_windows: int = 50,
    seed: int = 42
) -> Dict:
    """
    Run complete Layer 4 SHAPE analysis.
    
    Args:
        cache_dir: Path to data/_cache
        output_dir: Path to benchmarks/outputs/layer4_shape
        fse_start: FSE window start (1-based)
        fse_end: FSE window end (1-based, inclusive)
        core_start: FSE core start
        core_end: FSE core end
        n_genome_windows: Number of genome control windows
        seed: Random seed
    
    Returns:
        Summary dictionary
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    figures_dir = output_dir.parent.parent / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    
    print("=== Layer 4: SHAPE Agreement Analysis ===\n")
    
    # Load genome
    genome_path = cache_dir / "genome" / "NC_045512.2.fasta"
    with open(genome_path) as f:
        lines = f.readlines()
        genome = ''.join(line.strip() for line in lines[1:])
    
    # Extract FSE sequence
    fse_seq = genome[fse_start - 1 : fse_end]
    fse_seq_rna = fse_seq.replace('T', 'U')
    
    print(f"FSE sequence: {fse_start}-{fse_end} ({len(fse_seq)} nt)")
    print(f"Core: {core_start}-{core_end}\n")
    
    # Load SHAPE datasets
    print("Loading SHAPE datasets...")
    datasets = load_shape_datasets(cache_dir)
    
    for name, data in datasets.items():
        n_valid = np.sum(~np.isnan(data))
        print(f"  {name}: {n_valid}/{len(data)} positions")
    
    print()
    
    # Compute unpaired probabilities
    print("Computing unpaired probabilities...")
    unpaired_probs_window = compute_unpaired_probabilities(fse_seq_rna)
    
    core_len = core_end - core_start + 1
    unpaired_probs_core = unpaired_probs_window[:core_len]
    
    # Analyze each dataset
    results = {
        "fse_start": fse_start,
        "fse_end": fse_end,
        "core_start": core_start,
        "core_end": core_end,
        "window_length": len(fse_seq),
        "core_length": core_len,
        "datasets": {}
    }
    
    print("\n=== Per-Dataset Analysis ===\n")
    
    for dataset_name, full_reactivity in datasets.items():
        print(f"--- {dataset_name} ---")
        
        # Extract window and core
        window_react = extract_window(full_reactivity, fse_start, fse_end)
        core_react = extract_window(full_reactivity, core_start, core_end)
        
        # Normalize
        window_react_norm = normalize_reactivity(window_react, method="percentile_90")
        core_react_norm = normalize_reactivity(core_react, method="percentile_90")
        
        # Coverage
        window_coverage = np.sum(~np.isnan(window_react)) / len(window_react)
        core_coverage = np.sum(~np.isnan(core_react)) / len(core_react)
        
        print(f"  Coverage: window {window_coverage:.1%}, core {core_coverage:.1%}")
        
        # Spearman correlation - window
        valid_window = ~np.isnan(window_react_norm)
        if np.sum(valid_window) >= 10:
            rho_window, pval_window = spearmanr(
                unpaired_probs_window[valid_window],
                window_react_norm[valid_window]
            )
        else:
            rho_window, pval_window = np.nan, np.nan
        
        print(f"  Spearman (window): {rho_window:.3f} (p={pval_window:.2e})")
        
        # Spearman correlation - core
        valid_core = ~np.isnan(core_react_norm)
        if np.sum(valid_core) >= 10:
            rho_core, pval_core = spearmanr(
                unpaired_probs_core[valid_core],
                core_react_norm[valid_core]
            )
        else:
            rho_core, pval_core = np.nan, np.nan
        
        print(f"  Spearman (core): {rho_core:.3f} (p={pval_core:.2e})")
        
        # AUROC
        auroc_window = compute_auroc(unpaired_probs_window, window_react_norm)
        auroc_core = compute_auroc(unpaired_probs_core, core_react_norm)
        
        print(f"  AUROC (window): {auroc_window['auroc']:.3f} "
              f"({auroc_window['n_reactive']} reactive, {auroc_window['n_unreactive']} unreactive)")
        print(f"  AUROC (core): {auroc_core['auroc']:.3f}")
        
        # Store results
        results["datasets"][dataset_name] = {
            "window": {
                "spearman": rho_window,
                "spearman_pvalue": pval_window,
                "auroc": auroc_window["auroc"],
                "n_reactive": int(auroc_window["n_reactive"]),
                "n_unreactive": int(auroc_window["n_unreactive"]),
                "coverage": float(window_coverage),
                "n_valid": int(np.sum(valid_window))
            },
            "core": {
                "spearman": rho_core,
                "spearman_pvalue": pval_core,
                "auroc": auroc_core["auroc"],
                "n_reactive": int(auroc_core["n_reactive"]),
                "n_unreactive": int(auroc_core["n_unreactive"]),
                "coverage": float(core_coverage),
                "n_valid": int(np.sum(valid_core))
            }
        }
        
        print()
    
    # Dataset agreement
    print("=== Dataset Agreement ===\n")
    
    agreement = []
    dataset_names = list(datasets.keys())
    
    for i, name1 in enumerate(dataset_names):
        for name2 in dataset_names[i+1:]:
            react1 = extract_window(datasets[name1], fse_start, fse_end)
            react2 = extract_window(datasets[name2], fse_start, fse_end)
            
            react1_norm = normalize_reactivity(react1, method="percentile_90")
            react2_norm = normalize_reactivity(react2, method="percentile_90")
            
            valid = ~np.isnan(react1_norm) & ~np.isnan(react2_norm)
            
            if np.sum(valid) >= 10:
                rho, pval = spearmanr(react1_norm[valid], react2_norm[valid])
            else:
                rho, pval = np.nan, np.nan
            
            agreement.append({
                "dataset1": name1,
                "dataset2": name2,
                "spearman": rho,
                "pvalue": pval,
                "n_overlap": int(np.sum(valid))
            })
            
            print(f"  {name1} vs {name2}: {rho:.3f} (p={pval:.2e}, n={np.sum(valid)})")
    
    results["dataset_agreement"] = agreement
    
    # Tier analysis (for one representative in vitro dataset)
    print("\n=== Tier Reactivity Analysis (Incarnato in vitro) ===\n")
    
    window_react = extract_window(datasets["incarnato_invitro"], fse_start, fse_end)
    window_react_norm = normalize_reactivity(window_react, method="percentile_90")
    
    # Need pair probability matrix
    prob_matrix = compute_pair_probabilities(fse_seq_rna)
    
    tier_stats = analyze_tier_reactivity(unpaired_probs_window, prob_matrix, window_react_norm)
    
    for tier, stats in tier_stats.items():
        print(f"  {tier}: mean={stats['mean_reactivity']:.3f} ± {stats['std_reactivity']:.3f} (n={stats['n_positions']})")
    
    results["tier_analysis"] = tier_stats
    
    # SHAPE-directed folding (use Incarnato in vitro)
    print("\n=== SHAPE-Directed Folding ===\n")
    
    shape_result = shape_directed_folding(fse_seq_rna, window_react_norm)
    
    print(f"  MFE unconstrained: {shape_result['mfe_unconstrained']:.2f} kcal/mol")
    print(f"  MFE SHAPE: {shape_result['mfe_shape']:.2f} kcal/mol")
    print(f"  BP distance (MFE): {shape_result['bp_distance_mfe']}")
    print(f"  BP distance (MEA): {shape_result['bp_distance_mea']}")
    
    results["shape_directed_folding"] = {
        "dataset_used": "incarnato_invitro",
        "mfe_unconstrained": shape_result["mfe_unconstrained"],
        "mfe_shape": shape_result["mfe_shape"],
        "bp_distance_mfe": shape_result["bp_distance_mfe"],
        "bp_distance_mea": shape_result["bp_distance_mea"],
        "mfe_structure": shape_result["mfe_structure"],
        "mfe_structure_shape": shape_result["mfe_structure_shape"]
    }
    
    # Genome control
    print("\n=== Genome Control Analysis ===\n")
    print(f"Sampling {n_genome_windows} windows of length {len(fse_seq)}...")
    
    df_control = genome_control_analysis(
        genome, len(fse_seq), n_genome_windows, datasets, seed
    )
    
    # Compute percentiles for FSE
    fse_percentiles = {}
    
    for dataset_name in dataset_names:
        df_sub = df_control[df_control['dataset'] == dataset_name]
        fse_rho = results["datasets"][dataset_name]["window"]["spearman"]
        
        if not np.isnan(fse_rho):
            valid_rho = df_sub['spearman'].dropna()
            percentile = (valid_rho < fse_rho).sum() / len(valid_rho) * 100 if len(valid_rho) > 0 else np.nan
            fse_percentiles[dataset_name] = percentile
            
            print(f"  {dataset_name}: FSE at {percentile:.1f}th percentile")
        else:
            fse_percentiles[dataset_name] = np.nan
    
    results["genome_control"] = {
        "n_windows": n_genome_windows,
        "window_length": len(fse_seq),
        "fse_percentiles": fse_percentiles
    }
    
    # Save control data
    df_control.to_csv(output_dir / "genome_control.csv", index=False)
    
    # Generate figures
    print("\n=== Generating Figures ===\n")
    
    # SHAPE tracks
    position = np.arange(fse_start, fse_end + 1)
    datasets_norm_window = {}
    
    for name, full_react in datasets.items():
        window_react = extract_window(full_react, fse_start, fse_end)
        datasets_norm_window[name] = normalize_reactivity(window_react, method="percentile_90")
    
    plot_shape_tracks(
        position,
        unpaired_probs_window,
        datasets_norm_window,
        figures_dir / "layer4_shape_tracks.png",
        title=f"SARS-CoV-2 FSE (NC_045512.2:{fse_start}-{fse_end})"
    )
    print("  Saved: layer4_shape_tracks.png")
    
    # Genome control
    fse_spearman = {
        name: results["datasets"][name]["window"]["spearman"]
        for name in dataset_names
    }
    
    plot_genome_control(
        df_control,
        fse_spearman,
        figures_dir / "layer4_genome_control.png"
    )
    print("  Saved: layer4_genome_control.png")
    
    # Save results
    results_path = output_dir / "layer4_shape_results.json"
    
    with open(results_path, "w") as f:
        json.dump(results, f, indent=2)
    
    print(f"\n✓ Results saved to {results_path}")
    
    return results


def main():
    """Command-line interface."""
    import argparse
    
    parser = argparse.ArgumentParser(description="Layer 4: SHAPE Agreement Analysis")
    parser.add_argument("--cache-dir", type=Path, default=Path("data/_cache"),
                       help="Path to data cache")
    parser.add_argument("--output-dir", type=Path, default=Path("benchmarks/outputs/layer4_shape"),
                       help="Output directory")
    parser.add_argument("--n-windows", type=int, default=50,
                       help="Number of genome control windows")
    parser.add_argument("--seed", type=int, default=42,
                       help="Random seed")
    
    args = parser.parse_args()
    
    run_layer4_shape_analysis(
        args.cache_dir,
        args.output_dir,
        n_genome_windows=args.n_windows,
        seed=args.seed
    )


if __name__ == "__main__":
    main()
