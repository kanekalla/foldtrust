"""Comprehensive benchmark runner for all analysis layers."""

import json
import subprocess
import sys
from pathlib import Path
from typing import Dict

from foldtrust.benchmark.layer2_accuracy import run_layer2_benchmark
from foldtrust.benchmark.layer3_calibration import run_layer3_calibration
from foldtrust.benchmark.layer5_robustness import run_layer5_analysis
from foldtrust.benchmark.scoring import run_layer1_tests
from foldtrust.benchmark.shape import run_layer4_shape_analysis


def run_all_benchmarks(
    output_dir: Path,
    verbose: bool = True,
    use_full: bool = False,
) -> Dict:
    """
    Run complete benchmark analysis: layers 1, 2, 3, 4_shape, 5, then render tables.

    Args:
        output_dir: Output directory for all results
        verbose: Print progress messages
        use_full: Run full benchmark without subsampling

    Returns:
        Dictionary with results from all layers
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Create figures directory
    figures_dir = output_dir / "figures"
    figures_dir.mkdir(exist_ok=True)

    cache_dir = Path("data/_cache")
    cases_dir = Path("data/cases")

    results = {
        "output_dir": str(output_dir),
        "use_full": use_full,
        "layers": {},
    }

    # Layer 1: Scoring Correctness
    if verbose:
        print("\n" + "=" * 70)
        print("LAYER 1: SCORING CORRECTNESS")
        print("=" * 70)

    try:
        layer1_output = output_dir / "layer1"
        layer1_output.mkdir(parents=True, exist_ok=True)
        layer1 = run_layer1_tests(layer1_output)
        results["layers"]["layer1"] = layer1
        if verbose:
            print(
                f"✓ Layer 1 complete: {layer1.get('tests_passed', 0)}/{layer1.get('tests_total', 0)} tests passed"
            )
    except Exception as e:
        if verbose:
            print(f"✗ Layer 1 failed: {e}")
        results["layers"]["layer1"] = {"error": str(e)}

    # Layer 2: Structure Accuracy
    if verbose:
        print("\n" + "=" * 70)
        print("LAYER 2: STRUCTURE ACCURACY")
        print("=" * 70)

    try:
        layer2_output = output_dir / "layer2"
        layer2_output.mkdir(parents=True, exist_ok=True)
        layer2 = run_layer2_benchmark(
            cache_dir=cache_dir,
            output_dir=layer2_output,
            use_full=use_full,
            sample_size=None if use_full else 200,
        )
        results["layers"]["layer2"] = layer2
        if verbose:
            print(f"✓ Layer 2 complete: {layer2.get('n_structures', 0)} structures analyzed")
    except Exception as e:
        if verbose:
            print(f"✗ Layer 2 failed: {e}")
        results["layers"]["layer2"] = {"error": str(e)}

    # Layer 3: Calibration
    if verbose:
        print("\n" + "=" * 70)
        print("LAYER 3: CALIBRATION")
        print("=" * 70)

    try:
        layer3_output = output_dir / "layer3"
        layer3_output.mkdir(parents=True, exist_ok=True)
        layer3 = run_layer3_calibration(
            cache_dir=cache_dir,
            output_dir=layer3_output,
            use_full=use_full,
            sample_size=None if use_full else 200,
        )
        results["layers"]["layer3"] = layer3
        if verbose:
            print(
                f"✓ Layer 3 complete: ECE={layer3.get('ece', 0):.4f}, AUROC={layer3.get('auroc', 0):.4f}"
            )
    except Exception as e:
        if verbose:
            print(f"✗ Layer 3 failed: {e}")
        results["layers"]["layer3"] = {"error": str(e)}

    # Layer 4: SHAPE Agreement
    if verbose:
        print("\n" + "=" * 70)
        print("LAYER 4: SHAPE AGREEMENT")
        print("=" * 70)

    try:
        layer4_output = output_dir / "layer4_shape"
        layer4_output.mkdir(parents=True, exist_ok=True)
        layer4 = run_layer4_shape_analysis(cache_dir, layer4_output)
        results["layers"]["layer4_shape"] = {
            "fse_start": layer4.get("fse_start"),
            "fse_end": layer4.get("fse_end"),
            "n_datasets": len(layer4.get("datasets", {})),
        }
        if verbose:
            print(f"✓ Layer 4 complete: {len(layer4.get('datasets', {}))} datasets analyzed")
    except Exception as e:
        if verbose:
            print(f"✗ Layer 4 failed: {e}")
        results["layers"]["layer4_shape"] = {"error": str(e)}

    # Layer 5: Robustness
    if verbose:
        print("\n" + "=" * 70)
        print("LAYER 5: ROBUSTNESS")
        print("=" * 70)

    try:
        layer5_output = output_dir / "layer5"
        layer5_output.mkdir(parents=True, exist_ok=True)
        layer5 = run_layer5_analysis(cases_dir, layer5_output, cache_dir)
        results["layers"]["layer5"] = {"n_cases": len(layer5)}
        if verbose:
            print(f"✓ Layer 5 complete: {len(layer5)} cases analyzed")

        # Generate Layer 5 figures
        if verbose:
            print("  Generating Layer 5 figures...")
        from foldtrust.benchmark.layer5_figures import generate_all_layer5_figures

        generate_all_layer5_figures(layer5_output, figures_dir)
    except Exception as e:
        if verbose:
            print(f"✗ Layer 5 failed: {e}")
        results["layers"]["layer5"] = {"error": str(e)}

    # Render tables
    if verbose:
        print("\n" + "=" * 70)
        print("RENDERING TABLES")
        print("=" * 70)

    try:
        # Check if render scripts exist
        render_layer123_script = Path("scripts/render_layer123_tables.py")
        render_layer5_script = Path("scripts/render_layer5_tables.py")

        if render_layer123_script.exists():
            if verbose:
                print("  Running render_layer123_tables.py...")
            result = subprocess.run(
                [sys.executable, str(render_layer123_script), str(output_dir)],
                check=False,
                capture_output=True,
                text=True,
            )
            if result.returncode != 0:
                raise RuntimeError(f"render_layer123_tables.py failed: {result.stderr}")
            if verbose and result.stdout:
                print(result.stdout)

        if render_layer5_script.exists():
            if verbose:
                print("  Running render_layer5_tables.py...")
            result = subprocess.run(
                [sys.executable, str(render_layer5_script), str(output_dir)],
                check=False,
                capture_output=True,
                text=True,
            )
            if result.returncode != 0:
                raise RuntimeError(f"render_layer5_tables.py failed: {result.stderr}")
            if verbose and result.stdout:
                print(result.stdout)

        if verbose:
            print("✓ Tables rendered")
    except Exception as e:
        if verbose:
            print(f"✗ Table rendering failed: {e}")
        results["table_rendering"] = {"error": str(e)}

    # Save combined results
    summary_path = output_dir / "benchmark_summary.json"
    with open(summary_path, "w") as f:
        json.dump(results, f, indent=2)

    if verbose:
        print("\n" + "=" * 70)
        print("BENCHMARK COMPLETE")
        print("=" * 70)
        print(f"\nResults saved to: {output_dir}")
        print(f"Summary: {summary_path}")
        print(f"Figures: {figures_dir}")

    return results
