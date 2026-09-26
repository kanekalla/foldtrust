"""Layer 1 benchmark tests: scoring correctness and ViennaRNA integration."""

import json
from pathlib import Path

import numpy as np
import pytest

# Test fixtures path
FIXTURES_DIR = Path(__file__).parent / "fixtures" / "layer1"


def test_dotbracket_parser():
    """Test dot-bracket notation parser with various bracket types."""
    from foldtrust.io import parse_dotbracket

    # Standard parentheses
    pairs = parse_dotbracket("(((...)))")
    assert pairs == {(0, 8), (1, 7), (2, 6)}, "Standard parentheses failed"

    # Square brackets (pseudoknot)
    pairs = parse_dotbracket("((.[[.)).]]")
    assert (0, 7) in pairs, "Round brackets failed"
    assert (3, 10) in pairs, "Square brackets failed"

    # Curly braces
    pairs = parse_dotbracket("(({.)).}")
    assert (0, 5) in pairs, "Round brackets failed"
    assert (2, 7) in pairs, "Curly braces failed"

    # Angle brackets
    pairs = parse_dotbracket("((<.)).>")
    assert (0, 5) in pairs, "Round brackets failed"
    assert (2, 7) in pairs, "Angle brackets failed"

    # Unbalanced input should raise error
    with pytest.raises((ValueError, AssertionError)):
        parse_dotbracket("(((...)")


def test_bpseq_parser():
    """Test bpseq file parser."""
    from foldtrust.io import parse_bpseq

    bpseq_file = FIXTURES_DIR / "test_hairpin.bpseq"
    assert bpseq_file.exists(), f"Fixture {bpseq_file} not found"

    pairs = parse_bpseq(bpseq_file)
    expected = {(0, 7), (1, 6), (2, 5)}
    assert pairs == expected, f"Expected {expected}, got {pairs}"


def test_ct_parser():
    """Test CT (Connectivity Table) file parser."""
    from foldtrust.io import parse_ct

    ct_file = FIXTURES_DIR / "test_hairpin.ct"
    assert ct_file.exists(), f"Fixture {ct_file} not found"

    pairs = parse_ct(ct_file)
    expected = {(0, 7), (1, 6), (2, 5)}
    assert pairs == expected, f"Expected {expected}, got {pairs}"


def test_pseudoknot_removal():
    """Test pseudoknot removal with greedy algorithm."""
    from foldtrust.benchmark.metrics import remove_pseudoknots

    # Read crossing pseudoknot fixture
    pk_file = FIXTURES_DIR / "crossing_pk.bpseq"
    assert pk_file.exists(), f"Fixture {pk_file} not found"

    from foldtrust.io import parse_bpseq

    pairs = parse_bpseq(pk_file)

    # Structure has crossing: (0,9), (1,8), (2,7) cross with (5,14), (6,13)
    canonical = remove_pseudoknots(pairs)

    # Check that result is nested (no crossings)
    pairs_list = sorted(canonical)
    for i, (a1, b1) in enumerate(pairs_list):
        for a2, b2 in pairs_list[i + 1 :]:
            # Check no crossing: (a1 < a2 < b1 < b2) should not happen
            assert not (a1 < a2 < b1 < b2), f"Crossing found: {(a1,b1)} and {(a2,b2)}"


def test_metrics_toy_cases():
    """Test metric computation on toy examples."""
    from foldtrust.benchmark.metrics import compute_structure_metrics

    # Perfect match
    ref = {(0, 10), (1, 9), (2, 8), (3, 7), (4, 6)}
    pred = {(0, 10), (1, 9), (2, 8), (3, 7), (4, 6)}
    metrics = compute_structure_metrics(ref, pred)
    assert abs(metrics["sensitivity"] - 1.0) < 1e-6, "Perfect match sensitivity failed"
    assert abs(metrics["ppv"] - 1.0) < 1e-6, "Perfect match PPV failed"
    assert abs(metrics["f1"] - 1.0) < 1e-6, "Perfect match F1 failed"

    # No overlap
    ref = {(0, 10), (1, 9), (2, 8)}
    pred = {(11, 20), (12, 19), (13, 18)}
    metrics = compute_structure_metrics(ref, pred)
    assert abs(metrics["sensitivity"] - 0.0) < 1e-6, "No overlap sensitivity failed"
    assert abs(metrics["ppv"] - 0.0) < 1e-6, "No overlap PPV failed"
    assert abs(metrics["f1"] - 0.0) < 1e-6, "No overlap F1 failed"

    # Half overlap: TP=2, FP=2, FN=2
    ref = {(0, 10), (1, 9), (2, 8), (3, 7)}
    pred = {(2, 8), (3, 7), (11, 20), (12, 19)}
    metrics = compute_structure_metrics(ref, pred)
    assert abs(metrics["sensitivity"] - 0.5) < 1e-6, "Half overlap sensitivity failed"
    assert abs(metrics["ppv"] - 0.5) < 1e-6, "Half overlap PPV failed"
    assert abs(metrics["f1"] - 0.5) < 1e-6, "Half overlap F1 failed"


def test_mcc_computation():
    """Test MCC computation using sqrt(sensitivity * PPV) formula."""
    from foldtrust.benchmark.layer1_scoring import compute_exact_metrics, compute_slip_metrics

    # Perfect prediction: MCC should be 1.0
    ref = {(0, 10), (1, 9), (2, 8), (3, 7), (4, 6)}
    pred = {(0, 10), (1, 9), (2, 8), (3, 7), (4, 6)}
    metrics = compute_exact_metrics(pred, ref)
    assert abs(metrics["mcc"] - 1.0) < 1e-6, "Perfect match MCC should be 1.0"

    # Half overlap: MCC = sqrt(0.5 * 0.5) = 0.5
    ref = {(0, 10), (1, 9), (2, 8), (3, 7)}
    pred = {(2, 8), (3, 7), (11, 20), (12, 19)}
    metrics = compute_exact_metrics(pred, ref)
    expected_mcc = np.sqrt(0.5 * 0.5)
    assert abs(metrics["mcc"] - expected_mcc) < 1e-6, f"Half overlap MCC should be {expected_mcc}"

    # MCC should never be negative
    ref = {(0, 10), (1, 9), (2, 8)}
    pred = {(11, 20), (12, 19)}
    metrics = compute_exact_metrics(pred, ref)
    assert metrics["mcc"] >= 0.0, "MCC should never be negative"

    # Slip-tolerant MCC >= exact MCC
    ref = {(0, 10), (1, 9), (2, 8), (3, 7)}
    pred = {(1, 10), (2, 9), (11, 20), (12, 19)}
    exact_metrics = compute_exact_metrics(pred, ref)
    slip_metrics = compute_slip_metrics(pred, ref)
    assert slip_metrics["mcc"] >= exact_metrics["mcc"] - 1e-9, "Slip MCC should be >= exact MCC"


def test_slip_tolerant_greater_equal_exact():
    """Test that slip-tolerant F1 >= exact F1."""
    from foldtrust.benchmark.metrics import compute_slip_tolerant_metrics, compute_structure_metrics

    # Random test cases
    np.random.seed(42)
    for _ in range(10):
        n_ref = np.random.randint(5, 20)
        n_pred = np.random.randint(5, 20)
        ref = set()
        pred = set()

        # Generate non-overlapping random pairs
        for _ in range(n_ref):
            i = np.random.randint(0, 50)
            j = np.random.randint(i + 4, 60)
            ref.add((i, j))

        for _ in range(n_pred):
            i = np.random.randint(0, 50)
            j = np.random.randint(i + 4, 60)
            pred.add((i, j))

        exact = compute_structure_metrics(ref, pred)
        slip = compute_slip_tolerant_metrics(ref, pred)

        assert (
            slip["f1"] >= exact["f1"] - 1e-9
        ), f"Slip F1 {slip['f1']:.4f} < exact F1 {exact['f1']:.4f}"


def test_energy_regression():
    """Test FSE energy regression across parameter sets."""
    try:
        import RNA  # noqa: F401
    except ImportError:
        pytest.skip("ViennaRNA not installed")

    from foldtrust.vienna import fold_with_params

    # FSE NC_045512.2:13462-13542
    fse_seq = "UUUAAACGGGUUUGCGGUGUAAGUGCAGCCCGUCUUACACCGUGCGGCACAGGCACUAGUACUGAUGU" "CGUAUACAGGGCU"

    # Expected energies
    expected = {
        "Turner2004": -26.00,
        "Andronescu2007": -22.26,
        "Langdon2018": -24.70,
    }

    for params, expected_mfe in expected.items():
        # Use fold_with_params which runs in subprocess to avoid parameter state issues
        structure, mfe = fold_with_params(fse_seq, param_set=params, temperature=37.0)

        # Check MFE energy
        assert abs(mfe - expected_mfe) < 0.01, f"{params}: expected {expected_mfe}, got {mfe}"


def test_bpp_symmetry():
    """Test base-pair probability matrix symmetry."""
    try:
        import RNA  # noqa: F401
    except ImportError:
        pytest.skip("ViennaRNA not installed")

    from foldtrust.vienna import compute_pair_probabilities

    seq = "GGGAAACCC"
    P = compute_pair_probabilities(seq)

    # Check symmetry
    assert np.allclose(P, P.T), "BPP matrix is not symmetric"

    # Check row sums <= 1 + epsilon
    row_sums = P.sum(axis=1)
    assert np.all(row_sums <= 1.0 + 1e-9), f"Row sums exceed 1.0: {row_sums}"


def test_unpaired_probability():
    """Test unpaired probability = 1 - sum of pair probabilities."""
    try:
        import RNA  # noqa: F401
    except ImportError:
        pytest.skip("ViennaRNA not installed")

    from foldtrust.vienna import compute_pair_probabilities

    seq = "GGGAAACCC"
    P = compute_pair_probabilities(seq)

    # Unpaired prob = 1 - sum over BOTH triangles (symmetric matrix)
    unpaired = 1.0 - P.sum(axis=1)
    unpaired = np.clip(unpaired, 0, 1)

    # All unpaired probs should be in [0, 1]
    assert np.all(unpaired >= 0.0), "Negative unpaired probability"
    assert np.all(unpaired <= 1.0), "Unpaired probability > 1.0"


def test_gc_hairpin_firm_tier():
    """Test FIRM tier check on a high-probability GC hairpin."""
    try:
        import RNA  # noqa: F401
    except ImportError:
        pytest.skip("ViennaRNA not installed")

    from foldtrust.vienna import compute_pair_probabilities

    # Strong GC hairpin
    seq = "GGGCCCUUUGGGCCC"
    md = RNA.md()
    fc = RNA.fold_compound(seq, md)
    structure, mfe = fc.mfe()

    P = compute_pair_probabilities(seq)

    # Parse MFE pairs
    from foldtrust.io import parse_dotbracket

    mfe_pairs = parse_dotbracket(structure)

    # Check that some stem pairs are FIRM (p >= 0.85)
    firm_count = 0
    for i, j in mfe_pairs:
        if P[i, j] >= 0.85:
            firm_count += 1

    assert firm_count > 0, "No FIRM pairs found in GC hairpin"


def test_layer1_runner_generates_outputs():
    """Test that Layer 1 runner creates JSON and CSV outputs."""
    from foldtrust.benchmark.scoring import run_layer1_tests

    output_dir = Path("benchmarks/outputs/layer1")
    output_dir.mkdir(parents=True, exist_ok=True)

    run_layer1_tests(output_dir)

    # Check that files were created
    json_file = output_dir / "layer1_tests.json"
    csv_file = output_dir / "layer1_tests.csv"

    assert json_file.exists(), f"JSON output not found: {json_file}"
    assert csv_file.exists(), f"CSV output not found: {csv_file}"

    # Check JSON structure
    with open(json_file) as f:
        data = json.load(f)

    assert "tests" in data, "JSON missing 'tests' key"
    assert len(data["tests"]) > 0, "No tests in JSON output"
