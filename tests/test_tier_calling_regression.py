"""Regression tests for tier calling to prevent the reported bug."""

from foldtrust.vienna import compute_pair_probabilities, fold_mfe, parse_stems


def test_strong_gc_hairpin_yields_firm():
    """
    Test that a strong GC-rich hairpin yields FIRM tier labels.

    Regression test: validates that tier calling works correctly for stable stems
    with high base-pair probabilities.
    """
    # Strong GC-rich hairpin should have high pair probabilities
    sequence = "GCGCGCGCAAAAGCGCGCGC"

    structure, mfe_energy = fold_mfe(sequence)
    prob_matrix = compute_pair_probabilities(sequence)
    stems = parse_stems(structure, prob_matrix)

    # Should have at least one stem
    assert len(stems) > 0, "Expected at least one stem"

    # The main stem should be FIRM (high probability)
    main_stem = stems[0]
    assert main_stem["flag"] == "firm", f"Expected FIRM stem, got {main_stem['flag']}"
    assert main_stem["mean_prob"] > 0.85, f"Expected prob > 0.85, got {main_stem['mean_prob']:.4f}"

    # All pairs in the firm stem should have high individual probabilities
    for i, j in main_stem["pairs"]:
        pair_prob = prob_matrix[i, j]
        assert pair_prob > 0.8, f"Pair ({i+1}, {j+1}) has low prob {pair_prob:.4f} in FIRM stem"


def test_weak_hairpin_yields_floppy_or_soft():
    """Test that a weak hairpin with AU pairs yields SOFT or FLOPPY."""
    # Weak AU hairpin
    sequence = "AAAUUUAAAUUU"

    structure, mfe_energy = fold_mfe(sequence)
    prob_matrix = compute_pair_probabilities(sequence)
    stems = parse_stems(structure, prob_matrix)

    if len(stems) > 0:
        main_stem = stems[0]
        # Should not be FIRM for such a weak hairpin
        assert main_stem["flag"] in [
            "soft",
            "floppy",
        ], f"Expected SOFT or FLOPPY for weak AU hairpin, got {main_stem['flag']}"
        assert (
            main_stem["mean_prob"] < 0.85
        ), f"Expected prob < 0.85 for weak hairpin, got {main_stem['mean_prob']:.4f}"


def test_mixed_stability_structure():
    """Test a structure with both stable and unstable regions."""
    # GC-rich stem (stable) + AU-rich stem (less stable)
    sequence = "GCGCGCAAAAAAGCGCGCAAAAAUUUUAAAAAAAA"

    structure, mfe_energy = fold_mfe(sequence)
    prob_matrix = compute_pair_probabilities(sequence)
    stems = parse_stems(structure, prob_matrix)

    # Should have at least one stem
    assert len(stems) > 0

    # Check that tier labels correlate with mean probability
    for stem in stems:
        if stem["flag"] == "firm":
            assert stem["mean_prob"] >= 0.85
        elif stem["flag"] == "soft":
            assert 0.5 <= stem["mean_prob"] < 0.85
        elif stem["flag"] == "floppy":
            assert stem["mean_prob"] < 0.5


def test_sars2_fse_has_firm_stems():
    """
    Test that SARS-CoV-2 FSE has FIRM stems.

    The SARS-CoV-2 frameshifting element from data/cases/sars2-fse is a known
    strongly structured RNA with multiple stems including FIRM stems.
    Expected stems: soft 0.730, soft 0.847, soft 0.564, firm 0.946, firm 0.981
    """
    # SARS-CoV-2 FSE sequence (NC_045512.2:13462-13542)
    sequence = (
        "UUUAAACGGGUUUGCGGUGUAAGUGCAGCCCGUCUUACACCGUGCGGCACAGGCACUAGUACUGAUGUCGUAUACAGGGCUUUU"
    )

    structure, mfe_energy = fold_mfe(sequence)
    prob_matrix = compute_pair_probabilities(sequence)
    stems = parse_stems(structure, prob_matrix)

    # Should have multiple stems (FSE has 5 stems: 3 SOFT, 2 FIRM)
    assert len(stems) >= 4, f"Expected at least 4 stems, got {len(stems)}"

    # Should have at least two FIRM stems (the long pseudoknot stems)
    firm_stems = [s for s in stems if s["flag"] == "firm"]
    assert len(firm_stems) >= 2, f"Expected at least 2 FIRM stems, got {len(firm_stems)}"

    # The FIRM stems should have high probability (>= 0.85)
    for stem in firm_stems:
        assert (
            stem["mean_prob"] >= 0.85
        ), f"FIRM stem should have prob >= 0.85, got {stem['mean_prob']:.4f}"


def test_tier_thresholds_are_correct():
    """Test that tier threshold values match the documented thresholds."""
    # This tests the exact threshold values: FIRM >= 0.85, SOFT >= 0.5, FLOPPY < 0.5
    import numpy as np

    # Mock structures with known probabilities
    test_cases = [
        (0.95, "firm"),
        (0.85, "firm"),  # Boundary
        (0.84, "soft"),
        (0.65, "soft"),
        (0.50, "soft"),  # Boundary
        (0.49, "floppy"),
        (0.25, "floppy"),
        (0.01, "floppy"),
    ]

    sequence = "GCGCGCGCAAAAGCGCGCGC"
    structure, _ = fold_mfe(sequence)
    prob_matrix = compute_pair_probabilities(sequence)

    # Manually override probabilities to test threshold logic
    for target_prob, expected_flag in test_cases:
        # Create a mock probability matrix with uniform probability
        mock_matrix = np.full_like(prob_matrix, 0.01)

        # Set the main stem pairs to target probability
        stems_parsed = parse_stems(structure, mock_matrix)
        if len(stems_parsed) > 0:
            main_stem = stems_parsed[0]
            for i, j in main_stem["pairs"]:
                mock_matrix[i, j] = target_prob
                mock_matrix[j, i] = target_prob

        # Re-parse with the mock probabilities
        stems_with_mock = parse_stems(structure, mock_matrix)
        if len(stems_with_mock) > 0:
            result_flag = stems_with_mock[0]["flag"]
            assert (
                result_flag == expected_flag
            ), f"For prob={target_prob}, expected {expected_flag}, got {result_flag}"
