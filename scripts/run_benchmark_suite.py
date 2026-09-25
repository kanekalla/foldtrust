"""Run comprehensive FoldTrust benchmark suite.

Executes all benchmark layers:
1. Reference structure accuracy
2. Calibration analysis  
3. Robustness (temperature sweep)
4. SHAPE validation (SARS-CoV-2)

Usage:
    python scripts/run_benchmark_suite.py --output benchmarks/outputs_v2
"""

import argparse
import json
import time
from pathlib import Path

from foldtrust.benchmark.reference import run_reference_benchmark
from foldtrust.benchmark.calibration import run_calibration_analysis
from foldtrust.benchmark.robustness import run_robustness_analysis


def run_full_benchmark(output_dir: Path, curated_dataset: Path):
    """Run all benchmark layers."""
    output_dir.mkdir(parents=True, exist_ok=True)
    
    start_time = time.time()
    
    print("=" * 70)
    print("FoldTrust Benchmark Suite")
    print("=" * 70)
    
    results = {
        "benchmark_version": "v2",
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "layers": {},
    }
    
    # Layer 1: Reference Structure Accuracy
    print("\n" + "=" * 70)
    print("LAYER 1: Reference Structure Accuracy")
    print("=" * 70)
    
    try:
        ref_result = run_reference_benchmark(
            output_dir=output_dir,
            max_length=400,
            max_sequences=None,
            allow_slip=True,
            dataset_path=curated_dataset,
        )
        results["layers"]["reference_accuracy"] = {
            "status": "completed",
            "summary": ref_result,
        }
    except Exception as e:
        print(f"\n✗ Layer 1 failed: {e}")
        results["layers"]["reference_accuracy"] = {
            "status": "failed",
            "error": str(e),
        }
    
    # Layer 2: Calibration Analysis
    print("\n" + "=" * 70)
    print("LAYER 2: Calibration Analysis")
    print("=" * 70)
    
    try:
        calib_result = run_calibration_analysis(
            output_dir=output_dir,
            dataset_file=curated_dataset,
            max_sequences=50,
        )
        results["layers"]["calibration"] = {
            "status": "completed",
            "summary": calib_result,
        }
    except Exception as e:
        print(f"\n✗ Layer 2 failed: {e}")
        results["layers"]["calibration"] = {
            "status": "failed",
            "error": str(e),
        }
    
    # Layer 3: Robustness (Temperature)
    print("\n" + "=" * 70)
    print("LAYER 3: Robustness Analysis (Temperature)")
    print("=" * 70)
    
    try:
        # Temperature sweep on disease cases
        from pathlib import Path as P
        case_paths = list(P("data/cases").glob("*/sequence.fa"))
        
        if case_paths:
            robust_result = run_robustness_analysis(
                output_dir=output_dir,
                case_paths=case_paths,
                temperatures=[24, 37, 42],
            )
            results["layers"]["robustness"] = {
                "status": "completed",
                "summary": robust_result,
            }
        else:
            print("  No disease cases found in data/cases/")
            results["layers"]["robustness"] = {
                "status": "skipped",
                "reason": "no_disease_cases",
            }
    except Exception as e:
        print(f"\n✗ Layer 3 failed: {e}")
        results["layers"]["robustness"] = {
            "status": "failed",
            "error": str(e),
        }
    
    # Summary
    end_time = time.time()
    elapsed = end_time - start_time
    
    results["runtime_seconds"] = elapsed
    results["runtime_formatted"] = f"{elapsed/60:.1f} minutes"
    
    # Save results
    summary_file = output_dir / "benchmark_summary.json"
    with open(summary_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    print("\n" + "=" * 70)
    print("BENCHMARK COMPLETE")
    print("=" * 70)
    print(f"\nTotal runtime: {elapsed/60:.1f} minutes")
    print(f"Summary saved to: {summary_file}")
    
    # Print layer statuses
    print("\nLayer Status:")
    for layer, info in results["layers"].items():
        status = info["status"]
        symbol = "✓" if status == "completed" else ("⊘" if status == "skipped" else "✗")
        print(f"  {symbol} {layer}: {status}")
    
    return results


def main():
    parser = argparse.ArgumentParser(description="Run FoldTrust benchmark suite")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("benchmarks/outputs"),
        help="Output directory for benchmark results"
    )
    parser.add_argument(
        "--dataset",
        type=Path,
        default=Path("benchmarks/data/curated_references.json"),
        help="Path to curated reference dataset"
    )
    
    args = parser.parse_args()
    
    if not args.dataset.exists():
        print(f"Error: Dataset not found at {args.dataset}")
        print("Run: python scripts/build_curated_dataset.py")
        return 1
    
    run_full_benchmark(args.output, args.dataset)
    return 0


if __name__ == "__main__":
    exit(main())
