#!/usr/bin/env python3
"""Run FoldTrust Benchmark Layers 1-3.

Layer 1: Core validation tests
Layer 2: Accuracy on reference datasets (ArchiveII, bpRNA TS0, Rfam)
Layer 3: Calibration analysis
"""

import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from foldtrust.benchmark.layer1_scoring import run_layer1_tests
from foldtrust.benchmark.layer2_accuracy import run_layer2_benchmark
from foldtrust.benchmark.layer3_calibration import run_layer3_calibration


def main():
    cache_dir = Path("data/_cache")
    output_dir = Path("benchmarks/outputs")

    if not cache_dir.exists():
        print(f"Error: {cache_dir} not found. Run scripts/fetch_benchmark_data.py first.")
        return 1

    print("=" * 70)
    print("FoldTrust Benchmark: Layers 1-3")
    print("=" * 70)
    print()

    print("Layer 1: Core validation tests")
    print("-" * 70)
    layer1_output = output_dir / "layer1"
    layer1_results = run_layer1_tests(layer1_output)
    print(f"  ✓ Layer 1 complete. Results in {layer1_output}")
    print()

    print("Layer 2: Accuracy benchmark (ArchiveII, bpRNA TS0, Rfam)")
    print("-" * 70)
    layer2_output = output_dir / "layer2"
    layer2_df = run_layer2_benchmark(
        cache_dir=cache_dir, output_dir=layer2_output, max_length=500, sample_size=200
    )
    print(f"  ✓ Layer 2 complete. {len(layer2_df)} structures. Results in {layer2_output}")
    print()

    print("Layer 3: Calibration analysis")
    print("-" * 70)
    layer3_output = output_dir / "layer3"
    layer3_summary = run_layer3_calibration(
        cache_dir=cache_dir, layer2_output_dir=layer2_output, output_dir=layer3_output
    )
    print(f"  ✓ Layer 3 complete. ECE={layer3_summary['ece']:.4f}. Results in {layer3_output}")
    print()

    print("=" * 70)
    print("All layers complete!")
    print("=" * 70)
    return 0


if __name__ == "__main__":
    sys.exit(main())
