"""Regression tests for ft.tl ensemble computation.

Tests that ft.tl.compute_ensemble and ft.tl.compute_unpaired_probs
match the reference vienna.py implementation on the SARS-CoV-2 FSE.
"""

import numpy as np
import pytest

import foldtrust as ft
from foldtrust.vienna import compute_pair_probabilities, parse_stems


FSE_SEQ = "UUUAAACGGGUUUGCGGUGUAAGUGCAGCCCGUCUUACACCGUGCGGCACAGGCACUAGUACUGAUGUCGUAUACAGGGCU"

# Expected stem flags and mean probabilities from vienna.parse_stems on FSE MFE
EXPECTED_STEMS = [
    ("soft", 0.730),
    ("soft", 0.847),
    ("soft", 0.564),
    ("firm", 0.946),
    ("firm", 0.981),
]


def test_tl_pair_probs_equal_vienna():
    """Test that ft.tl.compute_ensemble pair probs equal vienna.compute_pair_probabilities."""
    fd = ft.FoldData(sequence=FSE_SEQ)
    ft.tl.compute_ensemble(fd)
    
    vienna_probs = compute_pair_probabilities(FSE_SEQ)
    
    assert fd.obsp['pair_probs'].shape == vienna_probs.shape
    np.testing.assert_allclose(
        fd.obsp['pair_probs'],
        vienna_probs,
        atol=1e-12,
        err_msg="ft.tl.compute_ensemble pair probs differ from vienna.compute_pair_probabilities"
    )


def test_tl_stem_flags_equal_vienna():
    """Test that ft.tl.call_tiers stem flags/means match vienna.parse_stems on FSE MFE."""
    fd = ft.FoldData(sequence=FSE_SEQ)
    
    try:
        ft.tl.fold_mfe(fd)
    except Exception:
        import RNA
        md = RNA.md()
        fc = RNA.fold_compound(FSE_SEQ, md)
        structure, energy = fc.mfe()
        fd.structures['mfe'] = structure
        fd.uns['mfe_energy'] = energy
    
    ft.tl.compute_ensemble(fd)
    ft.tl.call_tiers(fd)
    
    stems = fd.uns['stems']
    
    assert len(stems) == len(EXPECTED_STEMS), (
        f"Expected {len(EXPECTED_STEMS)} stems, got {len(stems)}"
    )
    
    for i, (stem, (expected_flag, expected_mean)) in enumerate(zip(stems, EXPECTED_STEMS)):
        assert stem['flag'] == expected_flag, (
            f"Stem {i}: expected flag '{expected_flag}', got '{stem['flag']}' "
            f"(mean_prob={stem['mean_prob']:.3f})"
        )
        np.testing.assert_allclose(
            stem['mean_prob'],
            expected_mean,
            atol=0.001,
            err_msg=f"Stem {i}: mean_prob differs from expected {expected_mean:.3f}"
        )


def test_unpaired_equals_one_minus_row_sums():
    """Test that ft.tl.compute_unpaired_probs equals 1 - row sums of symmetric pair_probs."""
    fd = ft.FoldData(sequence=FSE_SEQ)
    ft.tl.compute_ensemble(fd)
    ft.tl.compute_unpaired_probs(fd)
    
    prob_matrix = fd.obsp['pair_probs']
    unpaired = fd.obs['unpaired_prob']
    
    expected_unpaired = 1.0 - prob_matrix.sum(axis=1)
    expected_unpaired = np.clip(expected_unpaired, 0.0, 1.0)
    
    np.testing.assert_allclose(
        unpaired,
        expected_unpaired,
        atol=1e-12,
        err_msg="Unpaired probs should equal 1 - row sums"
    )
