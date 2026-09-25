"""Benchmark analysis for FoldTrust predictions."""

from foldtrust.benchmark.metrics import (
    compute_accuracy_metrics,
    compute_calibration_metrics,
    compute_probing_metrics,
)
from foldtrust.benchmark.reference import run_reference_benchmark
from foldtrust.benchmark.calibration import run_calibration_analysis
from foldtrust.benchmark.probing import run_probing_analysis
from foldtrust.benchmark.robustness import run_robustness_analysis

__all__ = [
    "compute_accuracy_metrics",
    "compute_calibration_metrics",
    "compute_probing_metrics",
    "run_reference_benchmark",
    "run_calibration_analysis",
    "run_probing_analysis",
    "run_robustness_analysis",
]
