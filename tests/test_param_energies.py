"""Test that ViennaRNA parameter sets produce different energies.

This is a regression test for the ViennaRNA parameter loading bug where
parameters appeared to not change because a fresh md() object was not created
after loading the parameter set.

Expected MFE energies for SARS-CoV-2 FSE (NC_045512.2:13462-13542):
- Turner2004:     -26.00 kcal/mol
- Andronescu2007: -22.26 kcal/mol
- Langdon2018:    -24.70 kcal/mol
"""

import pytest

try:
    import RNA
    HAS_RNA = True
except ImportError:
    HAS_RNA = False

from foldtrust.vienna import fold_with_params, ViennaRNAError


# SARS-CoV-2 FSE sequence (slippery site 13462-13468 through stem-loop 2)
FSE_SEQUENCE = (
    "UUUAAACGGGUUUGCGGUGUAAGUGCAGCCCGUCUUACACCGUGCGGCACAGGCACUAGUACUGAUG"
    "UCGUAUACAGGGCU"
)


@pytest.mark.skipif(not HAS_RNA, reason="ViennaRNA Python bindings not available")
def test_parameter_sets_differ():
    """Test that different parameter sets produce different MFE energies."""
    
    # Fold with each parameter set
    _, mfe_turner = fold_with_params(FSE_SEQUENCE, "Turner2004", 37.0)
    _, mfe_andronescu = fold_with_params(FSE_SEQUENCE, "Andronescu2007", 37.0)
    _, mfe_langdon = fold_with_params(FSE_SEQUENCE, "Langdon2018", 37.0)
    
    # Check expected values (allow 0.1 kcal/mol tolerance)
    assert abs(mfe_turner - (-26.00)) < 0.1, f"Turner2004: expected -26.00, got {mfe_turner:.2f}"
    assert abs(mfe_andronescu - (-22.26)) < 0.1, f"Andronescu2007: expected -22.26, got {mfe_andronescu:.2f}"
    assert abs(mfe_langdon - (-24.70)) < 0.1, f"Langdon2018: expected -24.70, got {mfe_langdon:.2f}"
    
    # Verify they are actually different
    assert mfe_turner != mfe_andronescu, "Turner2004 and Andronescu2007 should differ"
    assert mfe_turner != mfe_langdon, "Turner2004 and Langdon2018 should differ"
    assert mfe_andronescu != mfe_langdon, "Andronescu2007 and Langdon2018 should differ"
    
    print(f"✓ Turner2004:     {mfe_turner:.2f} kcal/mol")
    print(f"✓ Andronescu2007: {mfe_andronescu:.2f} kcal/mol")
    print(f"✓ Langdon2018:    {mfe_langdon:.2f} kcal/mol")


@pytest.mark.skipif(not HAS_RNA, reason="ViennaRNA Python bindings not available")
def test_invalid_parameter_set():
    """Test that invalid parameter set raises error."""
    
    with pytest.raises(ViennaRNAError, match="Unknown parameter set"):
        fold_with_params(FSE_SEQUENCE, "InvalidParams", 37.0)


@pytest.mark.skipif(not HAS_RNA, reason="ViennaRNA Python bindings not available")
def test_temperature_affects_energy():
    """Test that temperature changes affect MFE energy."""
    
    _, mfe_24 = fold_with_params(FSE_SEQUENCE, "Turner2004", 24.0)
    _, mfe_37 = fold_with_params(FSE_SEQUENCE, "Turner2004", 37.0)
    _, mfe_42 = fold_with_params(FSE_SEQUENCE, "Turner2004", 42.0)
    
    # Verify that energies differ with temperature
    # (exact values depend on thermodynamics, just check they're different)
    assert mfe_24 != mfe_37, "24°C and 37°C should produce different energies"
    assert mfe_37 != mfe_42, "37°C and 42°C should produce different energies"
    
    print(f"✓ MFE at 24°C: {mfe_24:.2f} kcal/mol")
    print(f"✓ MFE at 37°C: {mfe_37:.2f} kcal/mol")
    print(f"✓ MFE at 42°C: {mfe_42:.2f} kcal/mol")
