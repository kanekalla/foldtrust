"""Tests for benchmark metrics."""

import numpy as np

from foldtrust.benchmark.metrics import (
    compute_accuracy_metrics,
    compute_calibration_metrics,
    compute_tier_accuracy,
    parse_structure_pairs,
    stratify_by_length,
)


def test_parse_structure_pairs():
    """Test parsing base pairs from dot-bracket notation."""
    structure = "(((...)))"
    pairs = parse_structure_pairs(structure)
    assert (0, 8) in pairs
    assert (1, 7) in pairs
    assert (2, 6) in pairs
    assert len(pairs) == 3


def test_parse_structure_pairs_nested():
    """Test parsing nested structures."""
    structure = "((..((..))..((..))..))"
    pairs = parse_structure_pairs(structure)
    # Outer pair: (0,21), (1,20)
    # Inner left: (4,9), (5,8)
    # Inner right: (12,17), (13,16)
    assert len(pairs) == 6


def test_compute_accuracy_metrics_perfect():
    """Test accuracy metrics with perfect prediction."""
    predicted = "(((...)))"
    reference = "(((...)))"

    metrics = compute_accuracy_metrics(predicted, reference, allow_slip=False)

    assert metrics["sensitivity"] == 1.0
    assert metrics["ppv"] == 1.0
    assert metrics["f1"] == 1.0
    assert metrics["mcc"] == 1.0


def test_compute_accuracy_metrics_partial():
    """Test accuracy metrics with partial overlap."""
    predicted = "(((...)))"
    reference = "(((...)))"

    metrics = compute_accuracy_metrics(predicted, reference, allow_slip=False)

    assert metrics["tp"] == 3
    assert metrics["fp"] == 0
    assert metrics["fn"] == 0
    assert metrics["sensitivity"] == 1.0
    assert metrics["ppv"] == 1.0


def test_compute_accuracy_metrics_no_overlap():
    """Test accuracy metrics with no overlap."""
    predicted = "((..))...."
    reference = ".....((...))"

    metrics = compute_accuracy_metrics(predicted, reference, allow_slip=False)

    assert metrics["tp"] == 0
    assert metrics["sensitivity"] == 0.0
    assert metrics["ppv"] == 0.0
    assert metrics["f1"] == 0.0


def test_compute_calibration_metrics():
    """Test calibration metrics computation."""
    n = 10
    prob_matrix = np.random.rand(n, n)
    prob_matrix = (prob_matrix + prob_matrix.T) / 2
    np.fill_diagonal(prob_matrix, 0)

    reference_pairs = {(0, 9), (1, 8), (2, 7)}

    calib = compute_calibration_metrics(prob_matrix, reference_pairs, n_bins=5)

    assert "ece" in calib
    assert "auroc" in calib
    assert "auprc" in calib
    assert 0 <= calib["ece"] <= 1
    assert 0 <= calib["auroc"] <= 1
    assert 0 <= calib["auprc"] <= 1
    assert len(calib["bin_means"]) <= 5


def test_compute_tier_accuracy():
    """Test tier accuracy computation."""
    stems = [
        {"flag": "firm", "pairs": [(0, 10), (1, 9), (2, 8)]},
        {"flag": "soft", "pairs": [(3, 7), (4, 6)]},
        {"flag": "floppy", "pairs": [(11, 15)]},
    ]

    reference_pairs = {(0, 10), (1, 9), (2, 8), (3, 7)}

    tier_acc = compute_tier_accuracy(stems, reference_pairs)

    assert tier_acc["firm"]["ppv"] == 1.0
    assert tier_acc["firm"]["count"] == 3
    assert tier_acc["soft"]["ppv"] == 0.5
    assert tier_acc["soft"]["count"] == 2
    assert tier_acc["floppy"]["ppv"] == 0.0
    assert tier_acc["floppy"]["count"] == 1


def test_stratify_by_length():
    """Test length stratification."""
    sequences = ["A" * 25, "A" * 75, "A" * 150, "A" * 600]
    bins = [0, 50, 100, 200, 500]

    stratified = stratify_by_length(sequences, bins)

    assert "0-50" in stratified
    assert "50-100" in stratified
    assert "100-200" in stratified
    assert "500+" in stratified
    assert 0 in stratified["0-50"]
    assert 1 in stratified["50-100"]
    assert 2 in stratified["100-200"]
    assert 3 in stratified["500+"]
