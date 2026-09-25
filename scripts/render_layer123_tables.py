#!/usr/bin/env python3
"""
Render markdown tables from Layer 1-3 benchmark outputs.
Generates benchmarks/outputs/layer123_tables.md
"""

import json
from pathlib import Path
import pandas as pd


def render_layer2_tables(output_dir: Path) -> str:
    """Render Layer 2 accuracy tables."""
    layer2_dir = output_dir / "layer2"
    
    # Load data
    overall_df = pd.read_csv(layer2_dir / "layer2_summary_overall.csv")
    dataset_df = pd.read_csv(layer2_dir / "layer2_summary_dataset.csv")
    with open(layer2_dir / "layer2_summary.json") as f:
        summary = json.load(f)
    
    md = ["## Layer 2: Structure Accuracy\n"]
    
    # Overall metrics table (mean of per-structure F1)
    md.append("### Overall Metrics (Mean of Per-Structure)\n")
    md.append("| Method | Sensitivity | PPV | F1 | MCC |")
    md.append("|--------|-------------|-----|----|----|")
    
    for method in ['mfe', 'mea', 'centroid']:
        method_data = overall_df[overall_df['method'] == method]
        sen = method_data[method_data['metric'] == 'sensitivity']['mean'].values[0]
        ppv = method_data[method_data['metric'] == 'ppv']['mean'].values[0]
        f1 = method_data[method_data['metric'] == 'f1']['mean'].values[0]
        mcc = method_data[method_data['metric'] == 'mcc']['mean'].values[0]
        md.append(f"| {method.upper()} | {sen:.3f} | {ppv:.3f} | {f1:.3f} | {mcc:.3f} |")
    
    # Per-dataset table
    md.append("\n### Per-Dataset F1 (Mean of Per-Structure)\n")
    md.append("| Dataset | N | MFE | MEA | Centroid |")
    md.append("|---------|---|----|-----|---------")
    
    for dataset in ['ArchiveII', 'Rfam_seed', 'bpRNA_TS0']:
        subset = dataset_df[dataset_df['dataset'] == dataset]
        n = summary['datasets'][dataset]
        mfe_f1 = subset[subset['method'] == 'mfe']['f1_mean'].values[0]
        mea_f1 = subset[subset['method'] == 'mea']['f1_mean'].values[0]
        centroid_f1 = subset[subset['method'] == 'centroid']['f1_mean'].values[0]
        md.append(f"| {dataset} | {n} | {mfe_f1:.3f} | {mea_f1:.3f} | {centroid_f1:.3f} |")
    
    # Sample description
    md.append("\n### Sample Description\n")
    md.append(f"- **Total structures**: {summary['n_structures']}\n")
    md.append("- **Sampling**: Random seed 42, 200 structures per dataset (after ≤500 nt filter), re-seeded before each dataset\n")
    md.append("- **Pool sizes**: ArchiveII 3854, bpRNA TS0 1305, Rfam seed 468 (all ≤500 nt)\n")
    md.append("- **Duplicates**: Not removed (9 duplicate-sequence pairs exist in the 600)\n")
    md.append("- **Pseudoknot handling**: Greedy nested-only selection removes crossing pairs from reference\n")
    md.append("- **Matching**: Exact base-pair matching (i,j exact)\n")
    md.append("- **Slip tolerance**: Predicted pair (i,j) matches if any of (i±1,j), (i,j±1), or (i,j) is in reference; ")
    md.append("reference pair recovered if any of (i±1,j), (i,j±1), or (i,j) is predicted. ")
    md.append("Slip F1 ≥ exact F1 by construction.\n")
    md.append("- **Averaging**: F1 = mean of per-structure F1 (primary); bootstrap 95% CI with 1000 resamples\n")
    
    return "\n".join(md)


def render_layer3_tables(output_dir: Path) -> str:
    """Render Layer 3 calibration tables."""
    layer3_dir = output_dir / "layer3"
    
    # Load data
    with open(layer3_dir / "layer3_summary.json") as f:
        summary = json.load(f)
    
    cal_curve = pd.read_csv(layer3_dir / "layer3_calibration_curve.csv")
    tier_df = pd.read_csv(layer3_dir / "layer3_tier_summary.csv")
    
    md = ["\n## Layer 3: Calibration\n"]
    
    # Overall calibration metrics
    md.append("### Overall Calibration\n")
    md.append("| Metric | Value |")
    md.append("|--------|------:|")
    md.append(f"| ECE | {summary['ece']:.4f} |")
    
    # Compute ECE restricted to p>=0.5
    high_p_bins = cal_curve[cal_curve['bin_mean'] >= 0.5]
    if len(high_p_bins) > 0:
        total_high = high_p_bins['bin_count'].sum()
        ece_high = sum(
            (row['bin_count'] / total_high) * abs(row['bin_mean'] - row['bin_accuracy'])
            for _, row in high_p_bins.iterrows()
        )
    else:
        ece_high = 0.0
    md.append(f"| ECE (p≥0.5) | {ece_high:.4f} |")
    
    md.append(f"| AUROC | {summary['auroc']:.4f} |")
    md.append(f"| AUPRC | {summary['auprc']:.4f} |")
    md.append(f"| Candidates (p>{summary['prob_floor']}) | {summary['n_candidate_pairs']} |")
    md.append(f"| Positive pairs | {summary['n_positive_pairs']} ({100*summary['n_positive_pairs']/summary['n_candidate_pairs']:.2f}%) |")
    
    # Reliability per bin
    md.append("\n### Reliability Per Bin\n")
    md.append("| Bin | Mean Prob | Fraction Correct | Count |")
    md.append("|-----|-----------|------------------|-------|")
    
    bin_edges = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    for i, (_, row) in enumerate(cal_curve.iterrows()):
        lower = bin_edges[i]
        upper = bin_edges[i+1]
        bin_label = f"[{lower:.1f}, {upper:.1f}{')'  if upper < 1.0 else ']'}"
        md.append(f"| {bin_label} | {row['bin_mean']:.3f} | {row['bin_accuracy']:.3f} | {int(row['bin_count'])} |")
    
    # MFE tier PPVs
    md.append("\n### MFE-Pair Tier PPV\n")
    md.append("| Tier | Pooled PPV | N Structures | N Pairs |")
    md.append("|------|------------|--------------|---------|")
    
    for _, row in tier_df.iterrows():
        md.append(f"| {row['tier']} | {row['ppv_mean']:.3f} | {int(row['n_structures'])} | {int(row['total_pairs'])} |")
    
    md.append("\n**Tier definitions**: FIRM = p≥0.85, SOFT = 0.5≤p<0.85, FLOPPY = p<0.5. ")
    md.append("MFE pairs are binned by their own base-pair probability. ")
    md.append("Pooled PPV = (sum correct pairs) / (sum all pairs) across all structures in tier. ")
    md.append("N Structures = number of structures with ≥1 pair in that tier (note: same structure can contribute to multiple tiers).")
    
    return "\n".join(md)


def main():
    """Generate markdown tables from benchmark outputs."""
    output_dir = Path("benchmarks/outputs")
    
    # Check that required files exist
    required = [
        output_dir / "layer2" / "layer2_summary.json",
        output_dir / "layer2" / "layer2_summary_overall.csv",
        output_dir / "layer2" / "layer2_summary_dataset.csv",
        output_dir / "layer3" / "layer3_summary.json",
        output_dir / "layer3" / "layer3_calibration_curve.csv",
        output_dir / "layer3" / "layer3_tier_summary.csv",
    ]
    
    missing = [f for f in required if not f.exists()]
    if missing:
        print("ERROR: Missing required output files:")
        for f in missing:
            print(f"  - {f}")
        print("\nRun: python scripts/run_layers_1_2_3.py")
        return
    
    # Render tables
    md_lines = [
        "# Layer 1-3 Benchmark Tables",
        "",
        "Generated from benchmark outputs in `benchmarks/outputs/layer2/` and `layer3/`.",
        "",
    ]
    
    md_lines.append(render_layer2_tables(output_dir))
    md_lines.append(render_layer3_tables(output_dir))
    
    # Write output
    output_path = output_dir / "layer123_tables.md"
    with open(output_path, 'w') as f:
        f.write("\n".join(md_lines))
    
    print(f"✓ Tables written to {output_path}")


if __name__ == "__main__":
    main()
