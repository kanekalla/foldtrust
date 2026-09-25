"""
Layer 3: Calibration - The core FoldTrust claim about reliability tiers.

Analyzes pooled candidate pairs from Layer 2 sets: reliability diagrams,
ECE, AUROC, AUPRC, and tier analysis (FIRM/SOFT/FLOPPY PPV and coverage).
"""

import json
import gzip
from pathlib import Path
from typing import Dict, List, Set, Tuple, Optional
import numpy as np
import pandas as pd
from collections import defaultdict
from sklearn.metrics import roc_auc_score, average_precision_score, roc_curve, precision_recall_curve

from foldtrust.benchmark.layer2_accuracy import (
    load_archiveii,
    load_bprna_ts0,
    load_rfam_seed,
    remove_pseudoknots,
    bootstrap_ci
)
from foldtrust.benchmark.layer1_scoring import (
    parse_dotbracket,
    sequence_sha1,
)


def load_sample_ids(layer2_output_dir: Path) -> Dict[str, Set[str]]:
    """Load sample IDs from Layer 2's sample_ids.csv."""
    import csv
    
    sample_ids_path = layer2_output_dir / 'sample_ids.csv'
    
    if not sample_ids_path.exists():
        raise FileNotFoundError(
            f"sample_ids.csv not found at {sample_ids_path}. "
            "Run Layer 2 first to generate the sample."
        )
    
    # Build a dict mapping dataset to set of sha1s
    dataset_ids = defaultdict(set)
    
    with open(sample_ids_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            dataset_ids[row['dataset']].add(row['seq_sha1'])
    
    return dict(dataset_ids)


def filter_by_sample_ids(structures: List[Dict], dataset_name: str, sample_ids: Dict[str, Set[str]]) -> List[Dict]:
    """Filter structures to match Layer 2 sample IDs."""
    if dataset_name not in sample_ids:
        return structures
    
    allowed_sha1s = sample_ids[dataset_name]
    filtered = []
    
    for struct in structures:
        sha1 = sequence_sha1(struct['sequence'])
        if sha1 in allowed_sha1s:
            filtered.append(struct)
    
    return filtered


def compute_pair_probabilities(sequence: str, temperature: float = 37.0) -> np.ndarray:
    """
    Compute base-pair probability matrix with ViennaRNA.
    
    Returns:
        NxN matrix where [i, j] is probability of pairing (0-indexed)
    """
    import RNA
    
    # Load Turner2004 and create new md
    RNA.params_load_RNA_Turner2004()
    md = RNA.md()
    md.temperature = temperature
    
    fc = RNA.fold_compound(sequence, md)
    fc.pf()
    
    bpp = fc.bpp()
    n = len(sequence)
    prob_matrix = np.zeros((n, n))
    
    for i in range(1, n + 1):
        for j in range(i + 1, n + 1):
            if i < len(bpp) and j < len(bpp[i]):
                # ViennaRNA bpp is 1-indexed
                prob_matrix[i-1, j-1] = bpp[i][j]
                prob_matrix[j-1, i-1] = bpp[i][j]
    
    return prob_matrix


def collect_candidate_pairs(
    structures: List[Dict],
    prob_floor: float = 1e-3,
    temperature: float = 37.0
) -> Tuple[List[float], List[int]]:
    """
    Collect all candidate pairs (p > prob_floor) with their probabilities and labels.
    
    Returns:
        (y_prob, y_true) where y_true[i] = 1 if pair is in reference, 0 otherwise
    """
    y_prob = []
    y_true = []
    
    for i, struct in enumerate(structures):
        if (i + 1) % 50 == 0:
            print(f"  Processed {i + 1}/{len(structures)}...")
        
        sequence = struct['sequence']
        ref_pairs = struct['pairs']
        
        # Remove pseudoknots
        ref_nested, _ = remove_pseudoknots(ref_pairs)
        
        try:
            # Compute pair probabilities
            prob_matrix = compute_pair_probabilities(sequence, temperature)
            
            n = len(sequence)
            for i in range(n):
                for j in range(i + 1, n):
                    prob = prob_matrix[i, j]
                    if prob > prob_floor:
                        y_prob.append(prob)
                        y_true.append(1 if (i, j) in ref_nested else 0)
        
        except Exception as e:
            print(f"  Error processing {struct['name']}: {e}")
    
    return y_prob, y_true


def compute_calibration_curve(
    y_prob: np.ndarray,
    y_true: np.ndarray,
    n_bins: int = 10
) -> Tuple[List[float], List[float], List[int]]:
    """
    Compute calibration curve (reliability diagram data).
    
    Returns:
        (bin_means, bin_accuracies, bin_counts)
    """
    bin_edges = np.linspace(0, 1, n_bins + 1)
    bin_means = []
    bin_accs = []
    bin_counts = []
    
    for b in range(n_bins):
        if b == n_bins - 1:
            mask = (y_prob >= bin_edges[b]) & (y_prob <= bin_edges[b + 1])
        else:
            mask = (y_prob >= bin_edges[b]) & (y_prob < bin_edges[b + 1])
        
        if np.sum(mask) > 0:
            bin_mean = np.mean(y_prob[mask])
            bin_acc = np.mean(y_true[mask])
            bin_count = np.sum(mask)
            
            bin_means.append(bin_mean)
            bin_accs.append(bin_acc)
            bin_counts.append(int(bin_count))
    
    return bin_means, bin_accs, bin_counts


def compute_ece(
    y_prob: np.ndarray,
    y_true: np.ndarray,
    n_bins: int = 10
) -> float:
    """Compute Expected Calibration Error."""
    bin_edges = np.linspace(0, 1, n_bins + 1)
    ece = 0.0
    n_total = len(y_true)
    
    for b in range(n_bins):
        if b == n_bins - 1:
            mask = (y_prob >= bin_edges[b]) & (y_prob <= bin_edges[b + 1])
        else:
            mask = (y_prob >= bin_edges[b]) & (y_prob < bin_edges[b + 1])
        
        if np.sum(mask) > 0:
            bin_mean = np.mean(y_prob[mask])
            bin_acc = np.mean(y_true[mask])
            bin_count = np.sum(mask)
            
            ece += (bin_count / n_total) * abs(bin_mean - bin_acc)
    
    return ece


def analyze_tiers(
    structures: List[Dict],
    firm_threshold: float = 0.85,
    soft_threshold: float = 0.50,
    temperature: float = 37.0
) -> pd.DataFrame:
    """
    Analyze tier performance: PPV and coverage for FIRM/SOFT/FLOPPY pairs.
    
    Returns:
        DataFrame with tier statistics
    """
    tier_results = []
    
    for i, struct in enumerate(structures):
        if (i + 1) % 50 == 0:
            print(f"  Processed {i + 1}/{len(structures)}...")
        
        sequence = struct['sequence']
        ref_pairs = struct['pairs']
        dataset = struct['dataset']
        family = struct['family']
        
        # Remove pseudoknots
        ref_nested, _ = remove_pseudoknots(ref_pairs)
        
        try:
            # Compute pair probabilities
            prob_matrix = compute_pair_probabilities(sequence, temperature)
            
            # Get MFE structure
            import RNA
            RNA.params_load_RNA_Turner2004()
            md = RNA.md()
            md.temperature = temperature
            fc = RNA.fold_compound(sequence, md)
            mfe_structure, _ = fc.mfe()
            mfe_pairs = parse_dotbracket(mfe_structure)
            
            # Classify MFE pairs into tiers
            tier_stats = {
                'FIRM': {'tp': 0, 'total': 0},
                'SOFT': {'tp': 0, 'total': 0},
                'FLOPPY': {'tp': 0, 'total': 0}
            }
            
            for i, j in mfe_pairs:
                prob = prob_matrix[i, j]
                
                if prob >= firm_threshold:
                    tier = 'FIRM'
                elif prob >= soft_threshold:
                    tier = 'SOFT'
                else:
                    tier = 'FLOPPY'
                
                tier_stats[tier]['total'] += 1
                if (i, j) in ref_nested:
                    tier_stats[tier]['tp'] += 1
            
            # Compute tier PPVs
            for tier in ['FIRM', 'SOFT', 'FLOPPY']:
                total = tier_stats[tier]['total']
                tp = tier_stats[tier]['tp']
                ppv = tp / total if total > 0 else np.nan
                
                # Also compute what fraction of reference pairs this tier recovered
                coverage = tp / len(ref_nested) if len(ref_nested) > 0 else 0.0
                
                tier_results.append({
                    'name': struct['name'],
                    'dataset': dataset,
                    'family': family,
                    'tier': tier,
                    'ppv': ppv,
                    'coverage': coverage,
                    'n_pairs': total,
                    'n_correct': tp,
                    'n_ref_total': len(ref_nested)
                })
        
        except Exception as e:
            print(f"  Error processing {struct['name']}: {e}")
    
    return pd.DataFrame(tier_results)


def run_layer3_calibration(
    cache_dir: Path,
    output_dir: Path,
    max_length: int = 500,
    use_full: bool = False,
    sample_size: Optional[int] = None,
    prob_floor: float = 1e-3,
    n_bins: int = 10
) -> Dict:
    """
    Run Layer 3 calibration analysis.
    
    Args:
        cache_dir: Path to data/_cache
        output_dir: Output directory
        max_length: Max sequence length
        use_full: Use all sequences
        prob_floor: Minimum probability to include a pair
        n_bins: Number of bins for calibration curve
        
    Returns:
        Summary dict
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    
    if use_full:
        max_length = None
        print("Running FULL calibration (all sequences)")
    else:
        print(f"Running calibration with max_length={max_length}")
    
    # Load datasets
    print("\nLoading datasets...")
    archiveii = load_archiveii(cache_dir, max_length)
    bprna = load_bprna_ts0(cache_dir, max_length)
    rfam = load_rfam_seed(cache_dir, max_length)
    
    if sample_size:
        # Load sample IDs from Layer 2
        print("  Loading sample IDs from Layer 2...")
        layer2_output_dir = output_dir.parent / 'layer2'
        sample_ids = load_sample_ids(layer2_output_dir)
        
        # Filter to match Layer 2 sample
        archiveii = filter_by_sample_ids(archiveii, 'ArchiveII', sample_ids)
        bprna = filter_by_sample_ids(bprna, 'bpRNA_TS0', sample_ids)
        rfam = filter_by_sample_ids(rfam, 'Rfam_seed', sample_ids)
        
        print(f"  Loaded {len(archiveii)} ArchiveII, {len(bprna)} bpRNA, {len(rfam)} Rfam (matching Layer 2)")
    
    all_structures = archiveii + bprna + rfam
    print(f"Total structures: {len(all_structures)}")
    
    # Collect candidate pairs
    print(f"\nCollecting candidate pairs (p > {prob_floor})...")
    y_prob, y_true = collect_candidate_pairs(all_structures, prob_floor)
    
    y_prob = np.array(y_prob)
    y_true = np.array(y_true)
    
    print(f"  Total candidate pairs: {len(y_prob)}")
    print(f"  Positive pairs: {np.sum(y_true)} ({100 * np.mean(y_true):.2f}%)")
    
    # Compute calibration metrics
    print("\nComputing calibration metrics...")
    
    bin_means, bin_accs, bin_counts = compute_calibration_curve(y_prob, y_true, n_bins)
    ece = compute_ece(y_prob, y_true, n_bins)
    
    if np.sum(y_true) > 0:
        auroc = roc_auc_score(y_true, y_prob)
        auprc = average_precision_score(y_true, y_prob)
    else:
        auroc = 0.0
        auprc = 0.0
    
    print(f"  ECE: {ece:.4f}")
    print(f"  AUROC: {auroc:.4f}")
    print(f"  AUPRC: {auprc:.4f}")
    
    # Save calibration data
    calibration_df = pd.DataFrame({
        'bin_mean': bin_means,
        'bin_accuracy': bin_accs,
        'bin_count': bin_counts
    })
    calibration_df.to_csv(output_dir / 'layer3_calibration_curve.csv', index=False)
    
    # Tier analysis
    print("\nAnalyzing tier performance...")
    tier_df = analyze_tiers(all_structures)
    tier_df.to_csv(output_dir / 'layer3_tier_analysis.csv', index=False)
    
    # Compute tier summaries with bootstrap CIs
    print("\nComputing tier summaries with bootstrap CIs...")
    tier_summaries = []
    
    for tier in ['FIRM', 'SOFT', 'FLOPPY']:
        tier_subset = tier_df[tier_df['tier'] == tier]
        
        if len(tier_subset) > 0:
            # Remove NaN PPVs (structures with 0 pairs in this tier)
            ppv_values = tier_subset['ppv'].dropna().values
            cov_values = tier_subset['coverage'].dropna().values
            
            if len(ppv_values) > 0:
                ppv_mean, ppv_lower, ppv_upper = bootstrap_ci(ppv_values, seed=0)
                cov_mean, cov_lower, cov_upper = bootstrap_ci(cov_values, seed=0)
                
                tier_summaries.append({
                    'tier': tier,
                    'ppv_mean': ppv_mean,
                    'ppv_ci_lower': ppv_lower,
                    'ppv_ci_upper': ppv_upper,
                    'coverage_mean': cov_mean,
                    'coverage_ci_lower': cov_lower,
                    'coverage_ci_upper': cov_upper,
                    'n_structures': len(tier_subset),
                    'total_pairs': tier_subset['n_pairs'].sum()
                })
    
    tier_summary_df = pd.DataFrame(tier_summaries)
    tier_summary_df.to_csv(output_dir / 'layer3_tier_summary.csv', index=False)
    
    # Print tier summary
    print("\nTier Summary:")
    for _, row in tier_summary_df.iterrows():
        print(f"  {row['tier']:7s}: PPV = {row['ppv_mean']:.3f} "
              f"[{row['ppv_ci_lower']:.3f}, {row['ppv_ci_upper']:.3f}], "
              f"Coverage = {row['coverage_mean']:.3f}, "
              f"N = {row['n_structures']}, Pairs = {row['total_pairs']}")
    
    # Overall summary (convert numpy types to Python types for JSON)
    summary = {
        'n_structures': int(len(all_structures)),
        'n_candidate_pairs': int(len(y_prob)),
        'n_positive_pairs': int(np.sum(y_true)),
        'bootstrap_seed': 0,
        'bootstrap_iterations': 1000,
        'ece': float(ece),
        'auroc': float(auroc),
        'auprc': float(auprc),
        'n_bins': int(n_bins),
        'prob_floor': float(prob_floor),
        'tier_summary': [
            {k: float(v) if isinstance(v, (np.floating, np.integer)) else v 
             for k, v in tier.items()}
            for tier in tier_summaries
        ]
    }
    
    with open(output_dir / 'layer3_summary.json', 'w') as f:
        json.dump(summary, f, indent=2)
    
    print(f"\nResults saved to {output_dir}")
    
    return summary


if __name__ == '__main__':
    import sys
    
    cache_dir = Path('data/_cache')
    output_dir = Path('benchmarks/outputs/layer3')
    
    use_full = '--full' in sys.argv
    
    summary = run_layer3_calibration(cache_dir, output_dir, max_length=500, use_full=use_full)
    
    print("\n" + "=" * 60)
    print("LAYER 3 CALIBRATION SUMMARY")
    print("=" * 60)
    print(f"Candidate pairs: {summary['n_candidate_pairs']}")
    print(f"Positive pairs: {summary['n_positive_pairs']}")
    print(f"ECE: {summary['ece']:.4f}")
    print(f"AUROC: {summary['auroc']:.4f}")
    print(f"AUPRC: {summary['auprc']:.4f}")
