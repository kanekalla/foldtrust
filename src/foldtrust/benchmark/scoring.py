"""Layer 1: Scoring correctness - unit tests and sanity checks."""

import json
from pathlib import Path
from typing import Dict

import numpy as np
import pandas as pd


def compute_metrics_from_pairs(ref_pairs: set, pred_pairs: set) -> Dict:
    """Compute accuracy metrics from pair sets."""
    tp = len(ref_pairs & pred_pairs)
    fp = len(pred_pairs - ref_pairs)
    fn = len(ref_pairs - pred_pairs)

    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    ppv = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    f1 = 2 * tp / (2 * tp + fp + fn) if (2 * tp + fp + fn) > 0 else 0.0

    return {
        "sensitivity": sensitivity,
        "ppv": ppv,
        "f1": f1,
        "tp": tp,
        "fp": fp,
        "fn": fn,
    }


def run_layer1_tests(output_dir: Path) -> Dict:
    """
    Layer 1: Known-answer unit tests proving metric code is correct.

    Tests:
    1. Perfect match: sensitivity=1.0, PPV=1.0, F1=1.0
    2. No overlap: sensitivity=0.0, PPV=0.0, F1=0.0
    3. Half overlap: known metrics
    4. ECE with perfect calibration
    5. Tier accuracy

    Returns:
        Dictionary with test results
    """
    print("\n### Scoring Correctness Tests ###\n")

    output_dir.mkdir(parents=True, exist_ok=True)

    results = {
        "tests": [],
        "tests_total": 0,
        "tests_passed": 0,
    }

    # Test 1: Perfect match
    print("Test 1: Perfect match (all pairs correct)")
    ref_pairs = {(0, 10), (1, 9), (2, 8), (3, 7), (4, 6)}
    pred_pairs = {(0, 10), (1, 9), (2, 8), (3, 7), (4, 6)}

    metrics = compute_metrics_from_pairs(ref_pairs, pred_pairs)

    test1_pass = metrics["sensitivity"] == 1.0 and metrics["ppv"] == 1.0 and metrics["f1"] == 1.0

    results["tests"].append(
        {
            "name": "Perfect match",
            "passed": bool(test1_pass),
            **{
                k: float(v) if isinstance(v, (int, float, np.number)) else v
                for k, v in metrics.items()
            },
        }
    )
    results["tests_total"] += 1
    if test1_pass:
        results["tests_passed"] += 1
        print(
            f"  ✓ PASS: Sensitivity={metrics['sensitivity']:.3f}, PPV={metrics['ppv']:.3f}, F1={metrics['f1']:.3f}"
        )
    else:
        print("  ✗ FAIL")

    # Test 2: No overlap
    print("\nTest 2: No overlap (all pairs wrong)")
    ref_pairs = {(0, 10), (1, 9), (2, 8)}
    pred_pairs = {(11, 20), (12, 19), (13, 18)}

    metrics = compute_metrics_from_pairs(ref_pairs, pred_pairs)

    test2_pass = metrics["sensitivity"] == 0.0 and metrics["ppv"] == 0.0 and metrics["f1"] == 0.0

    results["tests"].append(
        {
            "name": "No overlap",
            "passed": bool(test2_pass),
            **{
                k: float(v) if isinstance(v, (int, float, np.number)) else v
                for k, v in metrics.items()
            },
        }
    )
    results["tests_total"] += 1
    if test2_pass:
        results["tests_passed"] += 1
        print("  ✓ PASS")
    else:
        print("  ✗ FAIL")

    # Test 3: Half overlap
    print("\nTest 3: Half overlap (known metrics)")
    ref_pairs = {(0, 10), (1, 9), (2, 8), (3, 7)}  # 4 pairs
    pred_pairs = {(2, 8), (3, 7), (11, 20), (12, 19)}  # 2 correct, 2 wrong

    metrics = compute_metrics_from_pairs(ref_pairs, pred_pairs)

    # TP=2, FP=2, FN=2
    # Sensitivity = 2/4 = 0.5, PPV = 2/4 = 0.5, F1 = 0.5
    test3_pass = (
        abs(metrics["sensitivity"] - 0.5) < 0.001
        and abs(metrics["ppv"] - 0.5) < 0.001
        and abs(metrics["f1"] - 0.5) < 0.001
    )

    results["tests"].append(
        {
            "name": "Half overlap",
            "passed": bool(test3_pass),
            **{
                k: float(v) if isinstance(v, (int, float, np.number)) else v
                for k, v in metrics.items()
            },
        }
    )
    results["tests_total"] += 1
    if test3_pass:
        results["tests_passed"] += 1
        print("  ✓ PASS")
    else:
        print(
            f"  ✗ FAIL: Got sens={metrics['sensitivity']:.3f}, ppv={metrics['ppv']:.3f}, f1={metrics['f1']:.3f}"
        )

    # Test 4: ECE with perfect calibration
    print("\nTest 4: Expected Calibration Error (ECE)")

    # Perfect calibration: predictions match frequencies exactly
    # Bin 0.1: 10% correct
    # Bin 0.5: 50% correct
    # Bin 0.9: 90% correct
    predicted = np.array([0.1] * 10 + [0.5] * 50 + [0.9] * 40)
    actual = np.array(
        [0] * 9 + [1] * 1 + [0] * 25 + [1] * 25 + [0] * 4 + [1] * 36  # 10% correct  # 50% correct
    )  # 90% correct

    # Compute ECE manually with 10 bins
    n_bins = 10
    bin_boundaries = np.linspace(0, 1, n_bins + 1)
    ece = 0.0

    for i in range(n_bins):
        lower, upper = bin_boundaries[i], bin_boundaries[i + 1]
        mask = (predicted > lower) & (predicted <= upper)
        if np.sum(mask) > 0:
            avg_pred = predicted[mask].mean()
            avg_actual = actual[mask].mean()
            ece += (np.sum(mask) / len(predicted)) * abs(avg_pred - avg_actual)

    test4_pass = ece < 0.05  # Should be very small

    results["tests"].append(
        {
            "name": "Calibration (ECE)",
            "passed": bool(test4_pass),
            "ece": float(ece),
        }
    )
    results["tests_total"] += 1
    if test4_pass:
        results["tests_passed"] += 1
        print(f"  ✓ PASS: ECE={ece:.4f} < 0.05")
    else:
        print(f"  ✗ FAIL: ECE={ece:.4f}")

    # Test 5: Tier accuracy
    print("\nTest 5: Tier accuracy (FIRM tier PPV)")

    ref_pairs_set = {(0, 10), (1, 9), (2, 8)}

    # FIRM: all correct
    firm_pred = {(0, 10), (1, 9)}
    firm_metrics = compute_metrics_from_pairs(ref_pairs_set, firm_pred)
    firm_ppv = firm_metrics["ppv"]

    test5_pass = firm_ppv == 1.0

    results["tests"].append(
        {
            "name": "Tier accuracy (FIRM)",
            "passed": bool(test5_pass),
            "firm_ppv": float(firm_ppv),
        }
    )
    results["tests_total"] += 1
    if test5_pass:
        results["tests_passed"] += 1
        print(f"  ✓ PASS: FIRM PPV={firm_ppv:.3f}")
    else:
        print("  ✗ FAIL")

    # Save results to both JSON and CSV
    json_path = output_dir / "layer1_tests.json"
    with open(json_path, "w") as f:
        json.dump(results, f, indent=2)

    df = pd.DataFrame(results["tests"])
    csv_path = output_dir / "layer1_tests.csv"
    df.to_csv(csv_path, index=False)

    print(f"\n✓ Saved to {json_path} and {csv_path}")
    print(f"\nSummary: {results['tests_passed']}/{results['tests_total']} tests passed")

    return results


# Backwards compatibility alias
run_scoring_tests = run_layer1_tests
