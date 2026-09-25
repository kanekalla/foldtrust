#!/usr/bin/env python3
"""
Run Layers 1-3 of the FoldTrust benchmark.

Layer 1: Scoring correctness tests
Layer 2: Structure accuracy on reference datasets
Layer 3: Calibration and tier analysis

Usage:
    python scripts/run_layers_1_2_3.py              # Default: max_length=500
    python scripts/run_layers_1_2_3.py --full       # All sequences (slow)
"""

import sys
import time
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.foldtrust.benchmark.layer1_scoring import run_layer1_tests
from src.foldtrust.benchmark.layer2_accuracy import run_layer2_benchmark
from src.foldtrust.benchmark.layer3_calibration import run_layer3_calibration
from src.foldtrust.benchmark.figures import generate_all_figures
import pandas as pd


def main():
    use_full = '--full' in sys.argv
    
    cache_dir = Path('data/_cache')
    base_output_dir = Path('benchmarks/outputs')
    
    if not cache_dir.exists():
        print("ERROR: data/_cache not found. Run scripts/fetch_benchmark_data.py first.")
        return 1
    
    print("=" * 70)
    print("FOLDTRUST BENCHMARK LAYERS 1-3")
    print("=" * 70)
    print()
    
    start_time = time.time()
    
    # Layer 1
    print("\n" + "=" * 70)
    print("LAYER 1: SCORING CORRECTNESS")
    print("=" * 70)
    
    layer1_output = base_output_dir / 'layer1'
    layer1_results = run_layer1_tests(layer1_output)
    
    # Layer 2
    print("\n" + "=" * 70)
    print("LAYER 2: STRUCTURE ACCURACY")
    print("=" * 70)
    
    layer2_output = base_output_dir / 'layer2'
    # Use smaller sample for default run (500 per dataset ~= 1500 total)
    # Full run with --full will use all ~5600 structures
    sample_size = None if use_full else 200
    layer2_summary = run_layer2_benchmark(
        cache_dir,
        layer2_output,
        max_length=500,
        use_full=use_full,
        sample_size=sample_size
    )
    
    # Layer 3
    print("\n" + "=" * 70)
    print("LAYER 3: CALIBRATION")
    print("=" * 70)
    
    layer3_output = base_output_dir / 'layer3'
    layer3_summary = run_layer3_calibration(
        cache_dir,
        layer3_output,
        max_length=500,
        use_full=use_full,
        sample_size=sample_size
    )
    
    # Generate figures
    print("\n" + "=" * 70)
    print("GENERATING FIGURES")
    print("=" * 70)
    
    figures_dir = base_output_dir / 'figures'
    figures_dir.mkdir(exist_ok=True)
    
    layer2_df = pd.read_csv(layer2_output / 'layer2_raw_results.csv')
    layer3_calib_df = pd.read_csv(layer3_output / 'layer3_calibration_curve.csv')
    layer3_tier_df = pd.read_csv(layer3_output / 'layer3_tier_summary.csv')
    
    generate_all_figures(layer2_df, layer3_calib_df, layer3_tier_df, figures_dir)
    
    # Final summary
    elapsed = time.time() - start_time
    
    print("\n" + "=" * 70)
    print("BENCHMARK COMPLETE")
    print("=" * 70)
    print(f"\nRuntime: {elapsed/60:.1f} minutes")
    print(f"\nResults saved to: {base_output_dir}")
    print(f"  Layer 1: {layer1_output}")
    print(f"  Layer 2: {layer2_output}")
    print(f"  Layer 3: {layer3_output}")
    print(f"  Figures: {figures_dir}")
    
    print("\n" + "=" * 70)
    print("HEADLINE METRICS")
    print("=" * 70)
    
    print(f"\nLayer 2 (Structure Accuracy):")
    print(f"  Total structures: {layer2_summary['n_structures']}")
    print(f"    ArchiveII: {layer2_summary['datasets']['ArchiveII']}")
    print(f"    bpRNA TS0: {layer2_summary['datasets']['bpRNA_TS0']}")
    print(f"    Rfam seed: {layer2_summary['datasets']['Rfam_seed']}")
    print(f"  MFE F1: {layer2_summary['mfe_f1_mean']:.3f} "
          f"[{layer2_summary['mfe_f1_ci'][0]:.3f}, {layer2_summary['mfe_f1_ci'][1]:.3f}]")
    print(f"  MEA F1: {layer2_summary['mea_f1_mean']:.3f} "
          f"[{layer2_summary['mea_f1_ci'][0]:.3f}, {layer2_summary['mea_f1_ci'][1]:.3f}]")
    
    print(f"\nLayer 3 (Calibration):")
    print(f"  Candidate pairs: {layer3_summary['n_candidate_pairs']:,}")
    print(f"  Positive pairs: {layer3_summary['n_positive_pairs']:,}")
    print(f"  ECE: {layer3_summary['ece']:.4f}")
    print(f"  AUROC: {layer3_summary['auroc']:.4f}")
    print(f"  AUPRC: {layer3_summary['auprc']:.4f}")
    
    print("\n  Tier PPVs (95% CI):")
    for tier_info in layer3_summary['tier_summary']:
        tier = tier_info['tier']
        ppv = tier_info['ppv_mean']
        ci_low = tier_info['ppv_ci_lower']
        ci_high = tier_info['ppv_ci_upper']
        print(f"    {tier:7s}: {ppv:.3f} [{ci_low:.3f}, {ci_high:.3f}]")
    
    print("\n" + "=" * 70)
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
