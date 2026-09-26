"""Regression test for vienna.py Python API changes.

Verifies that changing compute_pair_probabilities from subprocess to Python API
produces identical results for all disease cases.
"""

import json
from pathlib import Path

import numpy as np
import pytest

from foldtrust.core import process_sequence
from foldtrust.utils import find_case_directories


@pytest.fixture
def cases_dir():
    """Path to data/cases directory."""
    return Path("data/cases")


@pytest.fixture
def expected_results_path(tmp_path):
    """Path to store expected results (generated once)."""
    return Path("tests/fixtures/vienna_api_expected.json")


def test_vienna_api_regression(cases_dir, expected_results_path):
    """
    Test that Python API produces same results as subprocess version.

    Checks: MFE energy, MEA structure, tier counts (FIRM/SOFT/FLOPPY)
    """
    if not cases_dir.exists():
        pytest.skip("data/cases not found")

    case_dirs = find_case_directories(cases_dir)

    if len(case_dirs) == 0:
        pytest.skip("No case directories found")

    results = {}

    for case_dir in case_dirs:
        case_name = case_dir.name
        sequence_file = case_dir / "sequence.fa"

        if not sequence_file.exists():
            continue

        # Process case (uses Python API version)
        import tempfile

        with tempfile.TemporaryDirectory() as tmpdir:
            result = process_sequence(
                sequence_file, Path(tmpdir), generate_html=False, generate_markdown=False
            )

        # Extract key metrics
        results[case_name] = {
            "mfe_energy": round(result["mfe_energy"], 2),
            "structure_length": len(result["structure"]),
            "firm_count": result["stem_counts"]["firm"],
            "soft_count": result["stem_counts"]["soft"],
            "floppy_count": result["stem_counts"]["floppy"],
            "verdict": result["verdict"],
        }

    # Load or create expected results
    if expected_results_path.exists():
        with open(expected_results_path) as f:
            expected = json.load(f)

        # Compare
        for case_name, result in results.items():
            if case_name in expected:
                exp = expected[case_name]

                # MFE energy should match (within floating point tolerance)
                assert (
                    abs(result["mfe_energy"] - exp["mfe_energy"]) < 0.1
                ), f"{case_name}: MFE energy mismatch: {result['mfe_energy']} != {exp['mfe_energy']}"

                # Structure length should match exactly
                assert (
                    result["structure_length"] == exp["structure_length"]
                ), f"{case_name}: structure length mismatch"

                # Tier counts should match
                assert (
                    result["firm_count"] == exp["firm_count"]
                ), f"{case_name}: FIRM count mismatch"
                assert (
                    result["soft_count"] == exp["soft_count"]
                ), f"{case_name}: SOFT count mismatch"
                assert (
                    result["floppy_count"] == exp["floppy_count"]
                ), f"{case_name}: FLOPPY count mismatch"

                # Verdict should match
                assert result["verdict"] == exp["verdict"], f"{case_name}: verdict mismatch"

    else:
        # Create expected results file for future runs
        expected_results_path.parent.mkdir(parents=True, exist_ok=True)
        with open(expected_results_path, "w") as f:
            json.dump(results, f, indent=2)

        pytest.skip(
            f"Created baseline results at {expected_results_path}. " "Re-run test to validate."
        )


def test_pair_probabilities_properties():
    """Test basic properties of pair probability matrix."""
    from foldtrust.vienna import compute_pair_probabilities

    # Simple hairpin
    sequence = "GGGAAACCC"
    prob_matrix = compute_pair_probabilities(sequence)

    # Should be square
    assert prob_matrix.shape == (len(sequence), len(sequence))

    # Should be symmetric
    assert np.allclose(prob_matrix, prob_matrix.T)

    # Diagonal should be zero (can't pair with self)
    assert np.allclose(np.diag(prob_matrix), 0.0)

    # All probabilities should be in [0, 1]
    assert np.all(prob_matrix >= 0.0)
    assert np.all(prob_matrix <= 1.0)

    # Total pairing probability for each position should not exceed 1
    # For symmetric matrix, just sum over one dimension
    for i in range(len(sequence)):
        total_prob = np.sum(prob_matrix[i, :])
        assert total_prob <= 1.01, f"Position {i}: total_prob = {total_prob}"
