"""Tests for ViennaRNA integration."""

import pytest
import numpy as np

from foldtrust.vienna import check_viennarna, fold_mfe, compute_pair_probabilities, parse_stems


def test_check_viennarna():
    """Test ViennaRNA availability check."""
    result = check_viennarna()
    assert isinstance(result, bool)


@pytest.mark.skipif(not check_viennarna(), reason="ViennaRNA not installed")
def test_fold_mfe():
    """Test MFE folding with a simple hairpin."""
    sequence = "GGGAAACCC"
    structure, energy = fold_mfe(sequence)
    
    assert isinstance(structure, str)
    assert len(structure) == len(sequence)
    assert isinstance(energy, float)
    assert energy < 0
    assert "(" in structure
    assert ")" in structure


@pytest.mark.skipif(not check_viennarna(), reason="ViennaRNA not installed")
def test_compute_pair_probabilities():
    """Test pair probability computation."""
    sequence = "GGGAAACCC"
    prob_matrix = compute_pair_probabilities(sequence)
    
    assert prob_matrix.shape == (len(sequence), len(sequence))
    assert np.all(prob_matrix >= 0)
    assert np.all(prob_matrix <= 1)
    assert np.allclose(prob_matrix, prob_matrix.T)


@pytest.mark.skipif(not check_viennarna(), reason="ViennaRNA not installed")
def test_parse_stems():
    """Test stem parsing and classification."""
    sequence = "GGGAAACCC"
    structure = "(((...)))"
    
    prob_matrix = np.zeros((len(sequence), len(sequence)))
    prob_matrix[0, 8] = 0.9
    prob_matrix[8, 0] = 0.9
    prob_matrix[1, 7] = 0.9
    prob_matrix[7, 1] = 0.9
    prob_matrix[2, 6] = 0.9
    prob_matrix[6, 2] = 0.9
    
    stems = parse_stems(structure, prob_matrix, min_stem_length=2)
    
    assert len(stems) > 0
    assert all("flag" in s for s in stems)
    assert all(s["flag"] in ["firm", "soft", "floppy"] for s in stems)
    assert all("mean_prob" in s for s in stems)


def test_parse_stems_classification():
    """Test stem classification into firm/soft/floppy."""
    structure = "(((...)))"
    n = len(structure)
    
    prob_matrix_firm = np.zeros((n, n))
    prob_matrix_firm[0, 8] = 0.95
    prob_matrix_firm[8, 0] = 0.95
    prob_matrix_firm[1, 7] = 0.95
    prob_matrix_firm[7, 1] = 0.95
    prob_matrix_firm[2, 6] = 0.95
    prob_matrix_firm[6, 2] = 0.95
    
    stems = parse_stems(structure, prob_matrix_firm)
    if stems:
        assert stems[0]["flag"] == "firm"
    
    prob_matrix_soft = np.zeros((n, n))
    prob_matrix_soft[0, 8] = 0.7
    prob_matrix_soft[8, 0] = 0.7
    prob_matrix_soft[1, 7] = 0.7
    prob_matrix_soft[7, 1] = 0.7
    prob_matrix_soft[2, 6] = 0.7
    prob_matrix_soft[6, 2] = 0.7
    
    stems = parse_stems(structure, prob_matrix_soft)
    if stems:
        assert stems[0]["flag"] == "soft"
    
    prob_matrix_floppy = np.zeros((n, n))
    prob_matrix_floppy[0, 8] = 0.3
    prob_matrix_floppy[8, 0] = 0.3
    prob_matrix_floppy[1, 7] = 0.3
    prob_matrix_floppy[7, 1] = 0.3
    prob_matrix_floppy[2, 6] = 0.3
    prob_matrix_floppy[6, 2] = 0.3
    
    stems = parse_stems(structure, prob_matrix_floppy)
    if stems:
        assert stems[0]["flag"] == "floppy"
