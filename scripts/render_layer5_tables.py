"""Render Layer 5 summary tables from CSV outputs.

This script reads the Layer 5 CSV outputs and computes pooled retention metrics
(total retained / total reference stems across all 5 cases) as well as per-case means.
Output is written to benchmarks/outputs/layer5/layer5_tables.md for inclusion in docs.
"""

from pathlib import Path
import pandas as pd
import numpy as np


def compute_pooled_retention(df: pd.DataFrame, condition_col: str) -> pd.DataFrame:
    """
    Compute pooled retention: total retained / total reference stems.
    
    Args:
        df: DataFrame with columns [case, condition_col, tier, n_stems, retention]
        condition_col: Name of condition column ('temperature' or 'parameters' or 'flank_size')
        
    Returns:
        DataFrame with pooled retention per condition per tier
    """
    results = []
    
    for condition in df[condition_col].unique():
        cond_df = df[df[condition_col] == condition]
        
        for tier in ['FIRM', 'SOFT', 'FLOPPY']:
            tier_df = cond_df[cond_df['tier'] == tier]
            
            if len(tier_df) == 0:
                continue
            
            total_stems = tier_df['n_stems'].sum()
            total_retained = (tier_df['n_stems'] * tier_df['retention']).sum()
            pooled_retention = total_retained / total_stems if total_stems > 0 else None
            
            results.append({
                condition_col: condition,
                'tier': tier,
                'total_stems': int(total_stems),
                'pooled_retention': pooled_retention,
            })
    
    return pd.DataFrame(results)


def compute_per_case_mean(df: pd.DataFrame, condition_col: str) -> pd.DataFrame:
    """
    Compute per-case mean retention (average retention across cases).
    
    Args:
        df: DataFrame with columns [case, condition_col, tier, n_stems, retention]
        condition_col: Name of condition column
        
    Returns:
        DataFrame with per-case mean retention per condition per tier
    """
    results = []
    
    for condition in df[condition_col].unique():
        cond_df = df[df[condition_col] == condition]
        
        for tier in ['FIRM', 'SOFT', 'FLOPPY']:
            tier_df = cond_df[cond_df['tier'] == tier]
            
            if len(tier_df) == 0:
                continue
            
            mean_retention = tier_df['retention'].mean()
            n_cases = len(tier_df)
            
            results.append({
                condition_col: condition,
                'tier': tier,
                'n_cases': n_cases,
                'mean_retention': mean_retention,
            })
    
    return pd.DataFrame(results)


def render_tables(output_dir: Path) -> str:
    """
    Render all Layer 5 summary tables from CSVs.
    
    Returns:
        Markdown string with all tables
    """
    md = []
    md.append("# Layer 5: Computed Summary Tables from CSVs")
    md.append("")
    md.append("All numbers below are computed from the CSV outputs, not hand-typed.")
    md.append("")
    
    # Temperature sweep
    md.append("## Temperature Sweep (vs 37°C Turner2004 baseline)")
    md.append("")
    
    temp_df = pd.read_csv(output_dir / "temperature_stem_retention.csv")
    
    # Pooled retention
    pooled = compute_pooled_retention(temp_df, 'temperature')
    md.append("### Pooled Retention (total retained / total reference stems)")
    md.append("")
    
    pivot = pooled.pivot_table(index='tier', columns='temperature', values='pooled_retention')
    # Sort columns: baseline first, then temperatures
    col_order = ['37C_baseline'] + [c for c in pivot.columns if c != '37C_baseline']
    pivot = pivot[[c for c in col_order if c in pivot.columns]]
    
    md.append(pivot.to_markdown())
    md.append("")
    
    # Per-case mean
    per_case = compute_per_case_mean(temp_df, 'temperature')
    md.append("### Per-Case Mean Retention (average across 5 cases)")
    md.append("")
    
    pivot = per_case.pivot_table(index='tier', columns='temperature', values='mean_retention')
    pivot = pivot[[c for c in col_order if c in pivot.columns]]
    
    md.append(pivot.to_markdown())
    md.append("")
    
    # Parameter sweep
    md.append("## Parameter Set Sweep (vs Turner2004 baseline)")
    md.append("")
    
    params_df = pd.read_csv(output_dir / "parameters_stem_retention.csv")
    
    # Pooled retention
    pooled = compute_pooled_retention(params_df, 'parameters')
    md.append("### Pooled Retention")
    md.append("")
    
    pivot = pooled.pivot_table(index='tier', columns='parameters', values='pooled_retention')
    col_order = ['Turner2004_baseline'] + [c for c in pivot.columns if c != 'Turner2004_baseline']
    pivot = pivot[[c for c in col_order if c in pivot.columns]]
    
    md.append(pivot.to_markdown())
    md.append("")
    
    # Per-case mean
    per_case = compute_per_case_mean(params_df, 'parameters')
    md.append("### Per-Case Mean Retention")
    md.append("")
    
    pivot = per_case.pivot_table(index='tier', columns='parameters', values='mean_retention')
    pivot = pivot[[c for c in col_order if c in pivot.columns]]
    
    md.append(pivot.to_markdown())
    md.append("")
    
    # Window context
    md.append("## Window Context (vs 0-nt baseline)")
    md.append("")
    
    window_df = pd.read_csv(output_dir / "window_context_stem_retention.csv")
    
    # Pooled retention
    pooled = compute_pooled_retention(window_df, 'flank_size')
    md.append("### Pooled Retention")
    md.append("")
    
    pivot = pooled.pivot_table(index='tier', columns='flank_size', values='pooled_retention')
    
    md.append(pivot.to_markdown())
    md.append("")
    
    # Per-case mean
    per_case = compute_per_case_mean(window_df, 'flank_size')
    md.append("### Per-Case Mean Retention")
    md.append("")
    
    pivot = per_case.pivot_table(index='tier', columns='flank_size', values='mean_retention')
    
    md.append(pivot.to_markdown())
    md.append("")
    
    # MEA BP distances (window context)
    md.append("## Window Context: Mean MEA BP Distance")
    md.append("")
    
    metrics_df = pd.read_csv(output_dir / "window_context_metrics.csv")
    bp_dist_means = metrics_df.groupby('flank_size')['bp_distance_mea'].mean()
    
    md.append("| Flank Size (nt) | Mean BP Distance |")
    md.append("|-----------------|------------------|")
    for flank_size, mean_dist in bp_dist_means.items():
        md.append(f"| {flank_size} | {mean_dist:.1f} |")
    md.append("")
    
    return "\n".join(md)


def main():
    output_dir = Path("benchmarks/outputs/layer5")
    
    if not output_dir.exists():
        print(f"Error: {output_dir} does not exist")
        return 1
    
    # Render tables
    tables_md = render_tables(output_dir)
    
    # Save to file
    output_file = output_dir / "layer5_tables.md"
    with open(output_file, "w") as f:
        f.write(tables_md)
    
    print(f"✓ Rendered tables to {output_file}")
    
    # Also print to console
    print("\n" + tables_md)
    
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
