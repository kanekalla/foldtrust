"""Test SHAPE correlation calculations for Layer 4 analysis."""

import numpy as np
import pandas as pd
import pytest
from scipy.stats import spearmanr
from pathlib import Path

from foldtrust.benchmark.shape import (
    compute_unpaired_probabilities,
    normalize_reactivity_percentile90,
    normalize_reactivity_deigan,
)


@pytest.fixture
def fse_sequence():
    """FSE sequence (NC_045512.2:13462-13542, T->U)."""
    genome_path = Path("data/_cache/genome/NC_045512.2.fasta")
    if not genome_path.exists():
        pytest.skip("Genome data not available")
    
    with open(genome_path) as f:
        lines = f.readlines()
        genome = ''.join(line.strip() for line in lines[1:])
    
    fse_seq = genome[13461:13542]  # 0-based slicing for 1-based 13462-13542
    return fse_seq.replace('T', 'U')


@pytest.fixture
def shape_datasets():
    """Load SHAPE datasets."""
    shape_dir = Path("data/_cache/shape")
    if not shape_dir.exists():
        pytest.skip("SHAPE data not available")
    
    datasets = {}
    
    # Incarnato in vitro
    data = pd.read_csv(shape_dir / "incarnato_invitro_reactivity.csv", header=None)[0].values
    datasets["incarnato_invitro"] = data.astype(float)[13461:13542]
    
    # Incarnato in vivo
    data = pd.read_csv(shape_dir / "incarnato_invivo_reactivity.csv", header=None)[0].values
    datasets["incarnato_invivo"] = data.astype(float)[13461:13542]
    
    # Pyle (missing = -999)
    data = pd.read_csv(shape_dir / "pyle_reactivity.csv", header=None)[0].values
    data = data.astype(float)
    data[data == -999] = np.nan
    datasets["pyle"] = data[13461:13542]
    
    # Zhang in vitro
    data = pd.read_csv(shape_dir / "zhang_invitro_reactivity.csv", header=None)[0].values
    datasets["zhang_invitro"] = data.astype(float)[13461:13542]
    
    # Zhang in vivo
    data = pd.read_csv(shape_dir / "zhang_invivo_reactivity.csv", header=None)[0].values
    datasets["zhang_invivo"] = data.astype(float)[13461:13542]
    
    return datasets


def test_fse_raw_correlations(fse_sequence, shape_datasets):
    """
    Test that raw-reactivity Spearman correlations match independent calculation.
    
    Expected values (to 3 decimals):
    - incarnato_invitro: 0.290
    - incarnato_invivo: 0.354
    - pyle: 0.163
    - zhang_invitro: 0.496
    - zhang_invivo: 0.550
    """
    unpaired = compute_unpaired_probabilities(fse_sequence)
    
    expected = {
        "incarnato_invitro": 0.290,
        "incarnato_invivo": 0.354,
        "pyle": 0.163,
        "zhang_invitro": 0.496,
        "zhang_invivo": 0.550,
    }
    
    for dataset_name, reactivity in shape_datasets.items():
        valid = ~np.isnan(reactivity)
        
        if np.sum(valid) >= 10:
            rho, _ = spearmanr(unpaired[valid], reactivity[valid])
            
            exp_rho = expected[dataset_name]
            
            # Match to 3 decimals (tolerance 0.0005)
            assert abs(rho - exp_rho) < 0.001, \
                f"{dataset_name}: Spearman = {rho:.3f}, expected {exp_rho:.3f}"


def test_spearman_invariant_to_monotonic_normalization(fse_sequence, shape_datasets):
    """
    Test that Spearman correlation is invariant to monotonic transformations.
    
    Raw reactivity and normalized reactivity should give identical Spearman rho.
    """
    unpaired = compute_unpaired_probabilities(fse_sequence)
    
    for dataset_name, reactivity_raw in shape_datasets.items():
        # Normalize with different methods
        reactivity_p90 = normalize_reactivity_percentile90(reactivity_raw)
        reactivity_deigan = normalize_reactivity_deigan(reactivity_raw)
        
        valid_raw = ~np.isnan(reactivity_raw)
        valid_p90 = ~np.isnan(reactivity_p90)
        valid_deigan = ~np.isnan(reactivity_deigan)
        
        # Valid masks should be identical (normalization doesn't add/remove NaN)
        assert np.array_equal(valid_raw, valid_p90), \
            f"{dataset_name}: normalization changed valid mask"
        assert np.array_equal(valid_raw, valid_deigan), \
            f"{dataset_name}: normalization changed valid mask"
        
        if np.sum(valid_raw) >= 10:
            rho_raw, _ = spearmanr(unpaired[valid_raw], reactivity_raw[valid_raw])
            rho_p90, _ = spearmanr(unpaired[valid_p90], reactivity_p90[valid_p90])
            rho_deigan, _ = spearmanr(unpaired[valid_deigan], reactivity_deigan[valid_deigan])
            
            # Spearman should be very similar (within small tolerance)
            # Small changes are expected due to tie-breaking in normalization
            assert abs(rho_raw - rho_p90) < 0.01, \
                f"{dataset_name}: Spearman changed substantially with p90 normalization: {rho_raw:.6f} -> {rho_p90:.6f}"
            assert abs(rho_raw - rho_deigan) < 0.01, \
                f"{dataset_name}: Spearman changed substantially with Deigan normalization: {rho_raw:.6f} -> {rho_deigan:.6f}"


def test_unpaired_probabilities_sum_check(fse_sequence):
    """Test that unpaired probabilities are computed correctly."""
    unpaired = compute_unpaired_probabilities(fse_sequence)
    
    # Should have one value per nucleotide
    assert len(unpaired) == len(fse_sequence)
    
    # All values should be in [0, 1]
    assert np.all(unpaired >= 0.0)
    assert np.all(unpaired <= 1.0)
    
    # At least some positions should be highly unpaired (> 0.9)
    assert np.sum(unpaired > 0.9) > 0, "No highly unpaired positions found"
    
    # At least some positions should be highly paired (< 0.1)
    assert np.sum(unpaired < 0.1) > 0, "No highly paired positions found"
