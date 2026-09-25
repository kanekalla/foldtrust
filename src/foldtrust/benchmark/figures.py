"""Generate figures for benchmark layers."""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import Optional


# Set style
plt.style.use('seaborn-v0_8-darkgrid')
sns.set_palette("colorblind")


def plot_reliability_diagram(
    bin_means: list,
    bin_accs: list,
    bin_counts: list,
    output_path: Path,
    dataset_name: str = "All datasets"
):
    """Plot reliability diagram (calibration curve)."""
    fig, ax = plt.subplots(figsize=(8, 8))
    
    # Perfect calibration line
    ax.plot([0, 1], [0, 1], 'k--', label='Perfect calibration', linewidth=2)
    
    # Calibration curve
    ax.plot(bin_means, bin_accs, 'o-', markersize=10, linewidth=2, label='Observed')
    
    # Add bin counts as text
    for x, y, count in zip(bin_means, bin_accs, bin_counts):
        ax.annotate(f'{count}', (x, y), textcoords="offset points", 
                   xytext=(0, 10), ha='center', fontsize=8)
    
    ax.set_xlabel('Predicted pair probability', fontsize=12)
    ax.set_ylabel('Fraction of true pairs', fontsize=12)
    ax.set_title(f'Reliability Diagram - {dataset_name}', fontsize=14, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.05, 1.05)
    ax.set_aspect('equal')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"  Saved reliability diagram to {output_path}")


def plot_f1_by_family(
    df: pd.DataFrame,
    output_path: Path,
    method: str = 'mfe'
):
    """Plot F1 by family (box or violin plot)."""
    method_df = df[df['method'] == method].copy()
    
    # Filter to families with at least 5 structures
    family_counts = method_df['family'].value_counts()
    valid_families = family_counts[family_counts >= 5].index
    method_df = method_df[method_df['family'].isin(valid_families)]
    
    if len(method_df) == 0:
        print(f"  No families with >= 5 structures for {method}")
        return
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Sort families by median F1
    family_order = method_df.groupby('family')['f1'].median().sort_values(ascending=False).index
    
    sns.boxplot(data=method_df, x='family', y='f1', order=family_order, ax=ax)
    
    ax.set_xlabel('RNA Family', fontsize=12)
    ax.set_ylabel('F1 Score', fontsize=12)
    ax.set_title(f'F1 by Family ({method.upper()})', fontsize=14, fontweight='bold')
    ax.set_ylim(-0.05, 1.05)
    plt.xticks(rotation=45, ha='right')
    ax.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"  Saved F1 by family plot to {output_path}")


def plot_f1_vs_length(
    df: pd.DataFrame,
    output_path: Path,
    methods: list = ['mfe', 'mea', 'centroid']
):
    """Plot F1 vs sequence length."""
    fig, ax = plt.subplots(figsize=(10, 6))
    
    for method in methods:
        method_df = df[df['method'] == method]
        
        # Bin by length for smoothing
        length_bins = [0, 50, 100, 150, 200, 300, 500, 1000]
        binned = []
        
        for i in range(len(length_bins) - 1):
            mask = (method_df['length'] >= length_bins[i]) & (method_df['length'] < length_bins[i+1])
            if mask.sum() > 0:
                bin_center = (length_bins[i] + length_bins[i+1]) / 2
                mean_f1 = method_df[mask]['f1'].mean()
                std_f1 = method_df[mask]['f1'].std()
                binned.append((bin_center, mean_f1, std_f1))
        
        if binned:
            lengths, means, stds = zip(*binned)
            ax.errorbar(lengths, means, yerr=stds, marker='o', label=method.upper(), 
                       capsize=5, linewidth=2, markersize=8)
    
    ax.set_xlabel('Sequence Length (nt)', fontsize=12)
    ax.set_ylabel('F1 Score', fontsize=12)
    ax.set_title('F1 vs Sequence Length', fontsize=14, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(-0.05, 1.05)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"  Saved F1 vs length plot to {output_path}")


def plot_method_comparison(
    df: pd.DataFrame,
    output_path: Path,
    dataset: Optional[str] = None
):
    """Compare MFE vs MEA vs centroid."""
    if dataset:
        df = df[df['dataset'] == dataset].copy()
        title = f'Method Comparison - {dataset}'
    else:
        title = 'Method Comparison - All Datasets'
    
    fig, ax = plt.subplots(figsize=(8, 6))
    
    methods = ['mfe', 'mea', 'centroid']
    metrics = ['sensitivity', 'ppv', 'f1']
    
    x = np.arange(len(methods))
    width = 0.25
    
    for i, metric in enumerate(metrics):
        means = [df[df['method'] == m][metric].mean() for m in methods]
        ax.bar(x + i * width, means, width, label=metric.upper())
    
    ax.set_xlabel('Method', fontsize=12)
    ax.set_ylabel('Score', fontsize=12)
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.set_xticks(x + width)
    ax.set_xticklabels([m.upper() for m in methods])
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3, axis='y')
    ax.set_ylim(0, 1.0)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"  Saved method comparison to {output_path}")


def plot_tier_ppv(
    tier_summary_df: pd.DataFrame,
    output_path: Path
):
    """Plot tier PPV with confidence intervals."""
    fig, ax = plt.subplots(figsize=(8, 6))
    
    tiers = ['FIRM', 'SOFT', 'FLOPPY']
    x = np.arange(len(tiers))
    
    means = tier_summary_df['ppv_mean'].values
    ci_lower = tier_summary_df['ppv_ci_lower'].values
    ci_upper = tier_summary_df['ppv_ci_upper'].values
    
    errors = np.array([means - ci_lower, ci_upper - means])
    
    bars = ax.bar(x, means, color=['#2ecc71', '#f39c12', '#e74c3c'], alpha=0.7)
    ax.errorbar(x, means, yerr=errors, fmt='none', ecolor='black', capsize=10, linewidth=2)
    
    ax.set_xlabel('Reliability Tier', fontsize=12)
    ax.set_ylabel('PPV (Positive Predictive Value)', fontsize=12)
    ax.set_title('Tier PPV with 95% Bootstrap CIs', fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(tiers)
    ax.set_ylim(0, 1.0)
    ax.grid(True, alpha=0.3, axis='y')
    
    # Add value labels
    for i, (mean, lower, upper) in enumerate(zip(means, ci_lower, ci_upper)):
        ax.text(i, mean + 0.05, f'{mean:.3f}\n[{lower:.3f}, {upper:.3f}]', 
               ha='center', va='bottom', fontsize=9)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"  Saved tier PPV plot to {output_path}")


def generate_all_figures(
    layer2_results: pd.DataFrame,
    layer3_calibration: pd.DataFrame,
    layer3_tier_summary: pd.DataFrame,
    output_dir: Path
):
    """Generate all benchmark figures."""
    print("\nGenerating figures...")
    
    # Layer 2 figures
    plot_method_comparison(layer2_results, output_dir / 'layer2_method_comparison.png')
    plot_f1_vs_length(layer2_results, output_dir / 'layer2_f1_vs_length.png')
    
    # F1 by family (Rfam only)
    rfam_df = layer2_results[layer2_results['dataset'] == 'Rfam_seed']
    if len(rfam_df) > 0:
        plot_f1_by_family(rfam_df, output_dir / 'layer2_f1_by_family.png', method='mfe')
    
    # Layer 3 figures
    plot_reliability_diagram(
        layer3_calibration['bin_mean'].tolist(),
        layer3_calibration['bin_accuracy'].tolist(),
        layer3_calibration['bin_count'].tolist(),
        output_dir / 'layer3_reliability_diagram.png'
    )
    
    plot_tier_ppv(layer3_tier_summary, output_dir / 'layer3_tier_ppv.png')
    
    print("All figures generated!")
