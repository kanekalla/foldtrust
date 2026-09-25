"""
Layer 2: Structure Accuracy - Score ViennaRNA predictions against reference datasets.

Evaluates MFE, MEA, and centroid structures on ArchiveII, bpRNA TS0, Rfam seed,
and SPOT-RNA PDB TS1. Reports sensitivity/PPV/F1/MCC per dataset, family, and
length bin with bootstrap 95% CIs.
"""

import json
import gzip
from pathlib import Path
from typing import Dict, List, Set, Tuple, Optional
import numpy as np
import pandas as pd
from collections import defaultdict

from foldtrust.benchmark.layer1_scoring import (
    parse_dotbracket,
    parse_bpseq,
    compute_exact_metrics,
    compute_slip_metrics,
    fold_mfe,
    compute_bpp_matrix,
    remove_pseudoknots,
    sequence_sha1,
)


def save_sample_ids(archiveii: List[Dict], bprna: List[Dict], rfam: List[Dict], output_dir: Path):
    """Save sample IDs to CSV for Layer 3 consistency."""
    import csv
    
    sample_ids_path = output_dir / 'sample_ids.csv'
    
    with open(sample_ids_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['dataset', 'name', 'length', 'seq_sha1'])
        
        for struct in archiveii:
            writer.writerow(['ArchiveII', struct['name'], struct['length'], 
                           sequence_sha1(struct['sequence'])])
        for struct in bprna:
            writer.writerow(['bpRNA_TS0', struct['name'], struct['length'],
                           sequence_sha1(struct['sequence'])])
        for struct in rfam:
            writer.writerow(['Rfam_seed', struct['name'], struct['length'],
                           sequence_sha1(struct['sequence'])])
    
    print(f"  Saved sample IDs to {sample_ids_path}")


def load_archiveii(
    cache_dir: Path,
    max_length: Optional[int] = 500
) -> List[Dict]:
    """
    Load ArchiveII dataset from JSONL.
    
    Args:
        cache_dir: Path to data/_cache directory
        max_length: Maximum sequence length (None for all)
        
    Returns:
        List of dicts with name, sequence, pairs, family, length
    """
    jsonl_path = cache_dir / 'archiveII' / 'archiveII.jsonl.gz'
    
    if not jsonl_path.exists():
        raise FileNotFoundError(f"ArchiveII data not found at {jsonl_path}")
    
    structures = []
    
    with gzip.open(jsonl_path, 'rt') as f:
        for line in f:
            entry = json.loads(line)
            
            # Filter by length
            if max_length and entry['length'] > max_length:
                continue
            
            # Convert 1-indexed pairs to 0-indexed tuples
            pairs = set()
            for i, j in entry['pairs']:
                pairs.add((i-1, j-1) if i-1 < j-1 else (j-1, i-1))
            
            structures.append({
                'name': entry['name'],
                'sequence': entry['sequence'].replace('T', 'U'),
                'pairs': pairs,
                'family': entry['family'],
                'length': entry['length'],
                'dataset': 'ArchiveII'
            })
    
    return structures


def load_bprna_ts0(cache_dir: Path, max_length: Optional[int] = 500) -> List[Dict]:
    """Load bpRNA TS0 canonical pairs dataset."""
    ts0_dir = cache_dir / 'data' / 'bpRNA_dataset-canonicals' / 'TS0'
    
    if not ts0_dir.exists():
        raise FileNotFoundError(f"bpRNA TS0 not found at {ts0_dir}")
    
    structures = []
    
    for bpseq_file in sorted(ts0_dir.glob('*.bpseq')):
        sequence, structure = parse_bpseq(bpseq_file)
        pairs = parse_dotbracket(structure)
        
        # Filter by length
        if max_length and len(sequence) > max_length:
            continue
        
        structures.append({
            'name': bpseq_file.stem,
            'sequence': sequence,
            'pairs': pairs,
            'family': 'unknown',  # bpRNA doesn't have family labels
            'length': len(sequence),
            'dataset': 'bpRNA_TS0'
        })
    
    return structures


def load_rfam_seed(cache_dir: Path, max_length: Optional[int] = 500) -> List[Dict]:
    """Load Rfam seed consensus-derived structures."""
    jsonl_path = cache_dir / 'rfam' / 'rfam_seed_ss.jsonl'
    
    if not jsonl_path.exists():
        raise FileNotFoundError(f"Rfam seed data not found at {jsonl_path}")
    
    structures = []
    
    with open(jsonl_path) as f:
        for line in f:
            entry = json.loads(line)
            
            # Filter by length
            if max_length and len(entry['sequence']) > max_length:
                continue
            
            # Parse structure
            pairs = parse_dotbracket(entry['structure'])
            
            structures.append({
                'name': entry['seq_id'],
                'sequence': entry['sequence'],
                'pairs': pairs,
                'family': entry['family'],
                'length': len(entry['sequence']),
                'dataset': 'Rfam_seed'
            })
    
    return structures


def remove_pseudoknots(pairs: Set[Tuple[int, int]]) -> Tuple[Set[Tuple[int, int]], int]:
    """
    Remove pseudoknotted pairs, keeping only nested pairs.
    
    Returns:
        (nested_pairs, num_removed)
    """
    pairs_list = sorted(pairs)
    nested = []
    
    for pair_idx, pair_elem in enumerate(pairs_list):
        # Ensure we have a proper 2-tuple
        if not isinstance(pair_elem, tuple) or len(pair_elem) != 2:
            raise ValueError(f"Expected (i, j) tuple at index {pair_idx}, got {pair_elem} (type={type(pair_elem)})")
        a, b = pair_elem
        
        is_nested = True
        for j, (c, d) in enumerate(nested):
            # Check if (a, b) pseudoknots with (c, d)
            # Pseudoknot: c < a < d < b or a < c < b < d
            if (c < a < d < b) or (a < c < b < d):
                is_nested = False
                break
        
        if is_nested:
            nested.append((a, b))
    
    return set(nested), len(pairs) - len(nested)


def predict_structure_vienna(
    sequence: str,
    method: str = 'mfe',
    temperature: float = 37.0,
    gamma: float = 1.0
) -> Set[Tuple[int, int]]:
    """
    Predict structure with ViennaRNA.
    
    Args:
        sequence: RNA sequence
        method: 'mfe', 'mea', or 'centroid'
        temperature: Temperature in Celsius
        gamma: Gamma parameter for MEA (default 1.0)
        
    Returns:
        Set of predicted pairs (i, j) with i < j, 0-indexed
    """
    import RNA
    
    # Load Turner2004 parameters and create new md
    RNA.params_load_RNA_Turner2004()
    md = RNA.md()
    md.temperature = temperature
    
    fc = RNA.fold_compound(sequence, md)
    
    if method == 'mfe':
        structure, energy = fc.mfe()
    elif method == 'mea':
        # Need partition function for MEA
        fc.pf()
        structure, energy = fc.MEA(gamma)
    elif method == 'centroid':
        fc.pf()
        structure, energy = fc.centroid()
    else:
        raise ValueError(f"Unknown method: {method}")
    
    return parse_dotbracket(structure)


def bootstrap_ci(values: List[float], n_bootstrap: int = 1000, ci: float = 0.95, seed: int = 0) -> Tuple[float, float, float]:
    """
    Compute bootstrap confidence interval.
    
    Args:
        values: List of values to bootstrap
        n_bootstrap: Number of bootstrap iterations
        ci: Confidence level (default 0.95)
        seed: Random seed for reproducibility (default 0)
    
    Returns:
        (mean, lower_ci, upper_ci)
    """
    values = np.array(values)
    if len(values) == 0:
        return 0.0, 0.0, 0.0
    
    mean = np.mean(values)
    
    # Bootstrap resampling with seeded RNG
    rng = np.random.default_rng(seed)
    bootstrap_means = []
    for _ in range(n_bootstrap):
        sample = rng.choice(values, size=len(values), replace=True)
        bootstrap_means.append(np.mean(sample))
    
    alpha = 1 - ci
    lower = np.percentile(bootstrap_means, 100 * alpha / 2)
    upper = np.percentile(bootstrap_means, 100 * (1 - alpha / 2))
    
    return mean, lower, upper


def evaluate_dataset(
    structures: List[Dict],
    methods: List[str] = ['mfe', 'mea', 'centroid'],
    temperature: float = 37.0,
    gamma: float = 1.0
) -> pd.DataFrame:
    """
    Evaluate ViennaRNA predictions on a dataset.
    
    Returns:
        DataFrame with per-structure results
    """
    results = []
    
    for i, struct in enumerate(structures):
        if (i + 1) % 50 == 0:
            print(f"  Processed {i + 1}/{len(structures)}...")
        
        name = struct['name']
        sequence = struct['sequence']
        ref_pairs = struct['pairs']
        dataset = struct['dataset']
        family = struct['family']
        length = struct['length']
        
        # Remove pseudoknots from reference
        ref_nested, pk_removed = remove_pseudoknots(ref_pairs)
        
        for method in methods:
            try:
                pred_pairs = predict_structure_vienna(
                    sequence,
                    method=method,
                    temperature=temperature,
                    gamma=gamma
                )
                
                # Compute metrics (exact)
                metrics_exact = compute_exact_metrics(pred_pairs, ref_nested)
                
                # Compute metrics (slip-tolerant)
                metrics_slip = compute_slip_metrics(pred_pairs, ref_nested)
                
                results.append({
                    'name': name,
                    'dataset': dataset,
                    'family': family,
                    'length': length,
                    'method': method,
                    'pk_removed': pk_removed,
                    'ref_pairs': len(ref_nested),
                    'pred_pairs': len(pred_pairs),
                    # Exact metrics
                    'sensitivity': metrics_exact['sensitivity'],
                    'ppv': metrics_exact['ppv'],
                    'f1': metrics_exact['f1'],
                    'mcc': metrics_exact['mcc'],
                    # Slip-tolerant metrics
                    'sensitivity_slip': metrics_slip['sensitivity'],
                    'ppv_slip': metrics_slip['ppv'],
                    'f1_slip': metrics_slip['f1'],
                    'mcc_slip': metrics_slip['mcc'],
                })
                
            except Exception as e:
                print(f"  Error processing {name} with {method}: {e}")
    
    return pd.DataFrame(results)


def bin_by_length(df: pd.DataFrame, bins: List[int] = [0, 50, 100, 200, 500, 1000]) -> str:
    """Add length_bin column to DataFrame."""
    def get_bin(length):
        for i in range(len(bins) - 1):
            if bins[i] <= length < bins[i+1]:
                return f"{bins[i]}-{bins[i+1]}"
        return f"{bins[-1]}+"
    
    df['length_bin'] = df['length'].apply(get_bin)
    return df


def compute_summary_stats(df: pd.DataFrame, group_cols: List[str], metric: str = 'f1') -> pd.DataFrame:
    """Compute summary statistics with bootstrap CIs."""
    summaries = []
    
    for group_vals, group_df in df.groupby(group_cols):
        if not isinstance(group_vals, tuple):
            group_vals = (group_vals,)
        
        values = group_df[metric].values
        
        if len(values) > 0:
            mean, lower, upper = bootstrap_ci(values, n_bootstrap=1000, ci=0.95, seed=0)
            
            summary = dict(zip(group_cols, group_vals))
            summary.update({
                f'{metric}_mean': mean,
                f'{metric}_ci_lower': lower,
                f'{metric}_ci_upper': upper,
                'count': len(values)
            })
            summaries.append(summary)
    
    return pd.DataFrame(summaries)


def run_layer2_benchmark(
    cache_dir: Path,
    output_dir: Path,
    max_length: int = 500,
    use_full: bool = False,
    sample_size: Optional[int] = None
) -> Dict:
    """
    Run Layer 2 structure accuracy benchmark.
    
    Args:
        cache_dir: Path to data/_cache
        output_dir: Output directory for results
        max_length: Max sequence length (default 500)
        use_full: If True, use all sequences (ignores max_length)
        sample_size: If set, sample this many structures per dataset
        
    Returns:
        Summary dict
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    
    if use_full:
        max_length = None
        print("Running FULL benchmark (all sequences, no length cap)")
    else:
        print(f"Running benchmark with max_length={max_length}")
        if sample_size:
            print(f"  Sampling {sample_size} structures per dataset")
    
    # Load datasets
    print("\nLoading ArchiveII...")
    archiveii = load_archiveii(cache_dir, max_length)
    if sample_size and len(archiveii) > sample_size:
        import random
        random.seed(42)
        archiveii = random.sample(archiveii, sample_size)
    print(f"  Loaded {len(archiveii)} structures")
    
    print("\nLoading bpRNA TS0...")
    bprna = load_bprna_ts0(cache_dir, max_length)
    if sample_size and len(bprna) > sample_size:
        import random
        random.seed(42)
        bprna = random.sample(bprna, sample_size)
    print(f"  Loaded {len(bprna)} structures")
    
    print("\nLoading Rfam seed...")
    rfam = load_rfam_seed(cache_dir, max_length)
    if sample_size and len(rfam) > sample_size:
        import random
        random.seed(42)
        rfam = random.sample(rfam, sample_size)
    print(f"  Loaded {len(rfam)} structures")
    
    # Save sample IDs for Layer 3 consistency
    if sample_size:
        save_sample_ids(archiveii, bprna, rfam, output_dir)
    
    # Combine all structures
    all_structures = archiveii + bprna + rfam
    print(f"\nTotal structures: {len(all_structures)}")
    
    # Evaluate
    print("\nEvaluating ViennaRNA predictions...")
    print("(This may take several minutes)")
    results_df = evaluate_dataset(
        all_structures,
        methods=['mfe', 'mea', 'centroid'],
        temperature=37.0,
        gamma=1.0
    )
    
    # Save raw results
    results_csv = output_dir / 'layer2_raw_results.csv'
    results_df.to_csv(results_csv, index=False)
    print(f"\nRaw results saved to {results_csv}")
    
    # Add length bins
    results_df = bin_by_length(results_df)
    
    # Compute summaries
    print("\nComputing summary statistics...")
    
    # Per dataset and method
    summary_dataset = compute_summary_stats(
        results_df,
        ['dataset', 'method'],
        metric='f1'
    )
    summary_dataset.to_csv(output_dir / 'layer2_summary_dataset.csv', index=False)
    
    # Per family and method (for datasets with families)
    rfam_df = results_df[results_df['dataset'] == 'Rfam_seed']
    if len(rfam_df) > 0:
        summary_family = compute_summary_stats(
            rfam_df,
            ['family', 'method'],
            metric='f1'
        )
        summary_family.to_csv(output_dir / 'layer2_summary_family.csv', index=False)
    
    # Per length bin and method
    summary_length = compute_summary_stats(
        results_df,
        ['length_bin', 'method'],
        metric='f1'
    )
    summary_length.to_csv(output_dir / 'layer2_summary_length.csv', index=False)
    
    # Overall summary
    overall = []
    for method in ['mfe', 'mea', 'centroid']:
        method_df = results_df[results_df['method'] == method]
        for metric in ['sensitivity', 'ppv', 'f1', 'mcc']:
            mean, lower, upper = bootstrap_ci(method_df[metric].values, seed=0)
            overall.append({
                'method': method,
                'metric': metric,
                'mean': mean,
                'ci_lower': lower,
                'ci_upper': upper
            })
    
    overall_df = pd.DataFrame(overall)
    overall_df.to_csv(output_dir / 'layer2_summary_overall.csv', index=False)
    
    print(f"\nSummary statistics saved to {output_dir}")
    
    # Return headline numbers
    mfe_f1 = overall_df[(overall_df['method'] == 'mfe') & (overall_df['metric'] == 'f1')].iloc[0]
    mea_f1 = overall_df[(overall_df['method'] == 'mea') & (overall_df['metric'] == 'f1')].iloc[0]
    
    summary = {
        'n_structures': len(all_structures),
        'datasets': {
            'ArchiveII': len(archiveii),
            'bpRNA_TS0': len(bprna),
            'Rfam_seed': len(rfam)
        },
        'bootstrap_seed': 0,
        'bootstrap_iterations': 1000,
        'mfe_f1_mean': mfe_f1['mean'],
        'mfe_f1_ci': (mfe_f1['ci_lower'], mfe_f1['ci_upper']),
        'mea_f1_mean': mea_f1['mean'],
        'mea_f1_ci': (mea_f1['ci_lower'], mea_f1['ci_upper']),
    }
    
    # Save summary JSON
    with open(output_dir / 'layer2_summary.json', 'w') as f:
        json.dump(summary, f, indent=2)
    
    return summary


if __name__ == '__main__':
    import sys
    
    cache_dir = Path('data/_cache')
    output_dir = Path('benchmarks/outputs/layer2')
    
    use_full = '--full' in sys.argv
    
    summary = run_layer2_benchmark(cache_dir, output_dir, max_length=500, use_full=use_full)
    
    print("\n" + "=" * 60)
    print("LAYER 2 SUMMARY")
    print("=" * 60)
    print(f"Total structures: {summary['n_structures']}")
    print(f"  ArchiveII: {summary['datasets']['ArchiveII']}")
    print(f"  bpRNA TS0: {summary['datasets']['bpRNA_TS0']}")
    print(f"  Rfam seed: {summary['datasets']['Rfam_seed']}")
    print()
    print(f"MFE F1: {summary['mfe_f1_mean']:.3f} "
          f"(95% CI: [{summary['mfe_f1_ci'][0]:.3f}, {summary['mfe_f1_ci'][1]:.3f}])")
    print(f"MEA F1: {summary['mea_f1_mean']:.3f} "
          f"(95% CI: [{summary['mea_f1_ci'][0]:.3f}, {summary['mea_f1_ci'][1]:.3f}])")
