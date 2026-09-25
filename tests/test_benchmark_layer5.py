"""Tests for Layer 5 robustness analysis."""

import pytest
from pathlib import Path
import numpy as np

from foldtrust.benchmark.layer5_robustness import (
    CASE_COORDS,
    compute_mea_structure,
    compute_ensemble_defect,
    compute_basepair_distance,
    compute_stem_retention,
)
from foldtrust.vienna import fold_with_params, compute_pair_probs_with_params, parse_stems
from foldtrust.utils import read_fasta

try:
    import RNA
    HAS_RNA = True
except ImportError:
    HAS_RNA = False


@pytest.mark.skipif(not HAS_RNA, reason="ViennaRNA not available")
def test_case_coordinates_match_layer0():
    """Verify that case coordinates match Layer 0 verified values."""
    
    expected = {
        "sars2-fse": (81, "NC_045512.2", 13462, 13542),
        "smn2-iss-n1": (154, "NG_008728.1", 31999, 32152),
        "cftr-5utr": (200, "NM_000492.4", 1, 200),
        "mapt-e10": (183, "NG_007398.2", 120818, 121000),
        "hcv-ires-dii": (75, "AF009606.1", 44, 118),
    }
    
    for case_name, (length, accession, start, end) in expected.items():
        coords = CASE_COORDS[case_name]
        assert coords["length"] == length, f"{case_name}: length mismatch"
        assert coords["accession"] == accession, f"{case_name}: accession mismatch"
        assert coords["start"] == start, f"{case_name}: start mismatch"
        assert coords["end"] == end, f"{case_name}: end mismatch"


@pytest.mark.skipif(not HAS_RNA, reason="ViennaRNA not available")
def test_cases_load_from_files():
    """Verify cases can be loaded and match their declared coordinates."""
    
    cases_dir = Path("data/cases")
    if not cases_dir.exists():
        pytest.skip("data/cases not found")
    
    for case_name, coords in CASE_COORDS.items():
        seq_file = cases_dir / case_name / "sequence.fa"
        if not seq_file.exists():
            pytest.skip(f"{case_name} sequence.fa not found")
        
        _, sequence = read_fasta(seq_file)
        assert len(sequence) == coords["length"], \
            f"{case_name}: sequence length {len(sequence)} != {coords['length']}"


@pytest.mark.skipif(not HAS_RNA, reason="ViennaRNA not available")
def test_fse_37c_turner2004_energy():
    """FSE at 37°C Turner2004 should give -26.00 kcal/mol."""
    
    cases_dir = Path("data/cases")
    if not cases_dir.exists():
        pytest.skip("data/cases not found")
    
    seq_file = cases_dir / "sars2-fse" / "sequence.fa"
    if not seq_file.exists():
        pytest.skip("sars2-fse not found")
    
    _, sequence = read_fasta(seq_file)
    structure, mfe = fold_with_params(sequence, "Turner2004", 37.0)
    
    assert abs(mfe - (-26.00)) < 0.1, f"Expected -26.00, got {mfe:.2f}"


@pytest.mark.skipif(not HAS_RNA, reason="ViennaRNA not available")
def test_retention_of_self_is_one():
    """Retention of 37°C condition against itself should be 1.0."""
    
    # Simple test sequence
    sequence = "GGGAAACCC"
    
    # Compute reference
    ref_structure, _ = fold_with_params(sequence, "Turner2004", 37.0)
    ref_prob_matrix = compute_pair_probs_with_params(sequence, "Turner2004", 37.0)
    ref_stems = parse_stems(ref_structure, ref_prob_matrix)
    ref_mea_pairs, _ = compute_mea_structure(sequence, "Turner2004", 37.0)
    
    # Same condition
    cond_prob_matrix = compute_pair_probs_with_params(sequence, "Turner2004", 37.0)
    cond_mea_pairs, _ = compute_mea_structure(sequence, "Turner2004", 37.0)
    
    # Compute retention
    retention = compute_stem_retention(ref_stems, cond_prob_matrix, cond_mea_pairs)
    
    # All tiers with stems should have retention = 1.0
    for tier_name, tier_data in retention.items():
        if tier_data["count"] > 0:
            assert abs(tier_data["retention"] - 1.0) < 0.01, \
                f"{tier_name}: retention should be ~1.0, got {tier_data['retention']}"


@pytest.mark.skipif(not HAS_RNA, reason="ViennaRNA not available")
def test_mea_structure_is_paired():
    """MEA structure should produce valid pairs."""
    
    sequence = "GGGAAACCC"
    mea_pairs, mea_energy = compute_mea_structure(sequence, "Turner2004", 37.0)
    
    # Should have some pairs
    assert len(mea_pairs) > 0, "MEA should produce pairs"
    
    # All pairs should be valid (i < j)
    for i, j in mea_pairs:
        assert i < j, f"Invalid pair ({i}, {j})"
        assert 0 <= i < len(sequence)
        assert 0 <= j < len(sequence)


@pytest.mark.skipif(not HAS_RNA, reason="ViennaRNA not available")
def test_ensemble_defect_properties():
    """Ensemble defect should be non-negative."""
    
    sequence = "GGGAAACCC"
    structure, _ = fold_with_params(sequence, "Turner2004", 37.0)
    defect = compute_ensemble_defect(sequence, structure, "Turner2004", 37.0)
    
    assert defect >= 0, "Ensemble defect must be non-negative"
    assert defect <= len(sequence), "Ensemble defect cannot exceed sequence length"


@pytest.mark.skipif(not HAS_RNA, reason="ViennaRNA not available")
def test_basepair_distance_properties():
    """BP distance should be symmetric and non-negative."""
    
    pairs1 = [(0, 8), (1, 7), (2, 6)]
    pairs2 = [(0, 8), (1, 7)]
    
    dist = compute_basepair_distance(pairs1, pairs2)
    
    assert dist >= 0, "BP distance must be non-negative"
    assert dist == 1, "Should be 1 (one pair in pairs1 not in pairs2)"
    
    # Symmetry
    dist_reverse = compute_basepair_distance(pairs2, pairs1)
    assert dist == dist_reverse, "BP distance should be symmetric"


@pytest.mark.skipif(not HAS_RNA, reason="ViennaRNA not available")
def test_basepair_distance_identical_is_zero():
    """BP distance of identical pairs should be zero."""
    
    pairs = [(0, 8), (1, 7), (2, 6)]
    dist = compute_basepair_distance(pairs, pairs)
    
    assert dist == 0, "Distance of identical pairs should be 0"


def test_case_coords_complete():
    """All five cases should have coordinate definitions."""
    
    expected_cases = {"sars2-fse", "smn2-iss-n1", "cftr-5utr", "mapt-e10", "hcv-ires-dii"}
    assert set(CASE_COORDS.keys()) == expected_cases
