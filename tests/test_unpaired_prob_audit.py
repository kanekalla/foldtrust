"""Audit unpaired probability calculation against ViennaRNA."""

import numpy as np
import RNA

from foldtrust.benchmark.shape import compute_unpaired_probabilities
from foldtrust._core import FoldData
import foldtrust as ft


def test_unpaired_prob_vs_vienna_fse():
    """Test unpaired probability on FSE matches ViennaRNA directly."""
    # SARS-CoV-2 FSE sequence
    fse_seq = "UUUAAACGGGUUUGCGGUGUAAGUGCAGCCCGUCUUACACCGUGCGGCACAGGCACUAGUACUGAUGUCGUAUACAGGGCU"

    # Compute using our function
    our_unpaired = compute_unpaired_probabilities(fse_seq)

    # Compute using ViennaRNA directly
    fc = RNA.fold_compound(fse_seq)
    fc.pf()
    bpp = fc.bpp()
    n = len(fse_seq)

    # ViennaRNA's way: bpp is 1-indexed upper-triangular
    vienna_unpaired = np.zeros(n)
    for i in range(1, n + 1):
        paired_prob = 0.0
        # Sum over all j != i
        for j in range(1, n + 1):
            if i != j:
                # bpp[i][j] only valid for i < j in upper-triangular
                # For i > j, need to use bpp[j][i]
                if i < j:
                    paired_prob += bpp[i][j]
                else:
                    paired_prob += bpp[j][i]
        vienna_unpaired[i - 1] = 1.0 - paired_prob

    # Check they match within numerical tolerance
    assert np.allclose(our_unpaired, vienna_unpaired, atol=1e-6), (
        f"Unpaired probabilities don't match ViennaRNA!\n"
        f"Max difference: {np.max(np.abs(our_unpaired - vienna_unpaired))}\n"
        f"Our unpaired: {our_unpaired[:10]}\n"
        f"Vienna unpaired: {vienna_unpaired[:10]}"
    )


def test_unpaired_prob_vs_vienna_random():
    """Test unpaired probability on random sequence matches ViennaRNA."""
    # Random 120-nt sequence
    np.random.seed(42)
    bases = ["A", "C", "G", "U"]
    seq = "".join(np.random.choice(bases, 120))

    # Compute using our function
    our_unpaired = compute_unpaired_probabilities(seq)

    # Compute using ViennaRNA directly
    fc = RNA.fold_compound(seq)
    fc.pf()
    bpp = fc.bpp()
    n = len(seq)

    vienna_unpaired = np.zeros(n)
    for i in range(1, n + 1):
        paired_prob = 0.0
        for j in range(1, n + 1):
            if i != j:
                if i < j:
                    paired_prob += bpp[i][j]
                else:
                    paired_prob += bpp[j][i]
        vienna_unpaired[i - 1] = 1.0 - paired_prob

    assert np.allclose(our_unpaired, vienna_unpaired, atol=1e-6), (
        f"Unpaired probabilities don't match ViennaRNA on random sequence!\n"
        f"Max difference: {np.max(np.abs(our_unpaired - vienna_unpaired))}"
    )


def test_unpaired_prob_symmetry():
    """Test that bpp access is symmetric: bpp[min(i,j)][max(i,j)] works."""
    seq = "GCGCAAAAGCGC"  # Simple hairpin

    fc = RNA.fold_compound(seq)
    fc.pf()
    bpp = fc.bpp()
    n = len(seq)

    # Verify symmetry: bpp[i][j] == bpp[j][i] when accessed correctly
    for i in range(1, n + 1):
        for j in range(i + 1, n + 1):  # Only check upper triangle
            # Our access method: min/max
            val1 = bpp[min(i, j)][max(i, j)]
            # Direct access (should be same since i < j)
            val2 = bpp[i][j]
            assert abs(val1 - val2) < 1e-10, f"bpp[{i}][{j}] access mismatch"

            # Reverse access should also work
            val3 = bpp[min(j, i)][max(j, i)]
            assert abs(val1 - val3) < 1e-10, f"bpp[{j}][{i}] reverse access mismatch"


def test_unpaired_prob_sum_constraint():
    """Test that unpaired_prob + sum(pair_probs) ≈ 1 for each position."""
    seq = "UUUAAACGGGUUUGCGGUGUAAGUGCAGCCCGUCUUACACCGUGCGGCACAGGCACUAGUACUGAUGUCGUAUACAGGGCU"

    fc = RNA.fold_compound(seq)
    fc.pf()
    bpp = fc.bpp()
    n = len(seq)

    unpaired = compute_unpaired_probabilities(seq)

    for i in range(1, n + 1):
        # Sum all pairing probabilities for position i
        paired_prob = sum(bpp[min(i, j)][max(i, j)] for j in range(1, n + 1) if i != j)
        total = unpaired[i - 1] + paired_prob

        # Should sum to 1 (within numerical precision)
        assert abs(total - 1.0) < 1e-6, f"Position {i}: unpaired + paired = {total:.10f} != 1.0"


def test_unpaired_prob_tl_path_fse():
    """Test ft.tl.compute_unpaired_probs path on FSE matches ViennaRNA."""
    # SARS-CoV-2 FSE sequence
    fse_seq = "UUUAAACGGGUUUGCGGUGUAAGUGCAGCCCGUCUUACACCGUGCGGCACAGGCACUAGUACUGAUGUCGUAUACAGGGCU"

    # Create FoldData and compute ensemble
    fd = FoldData(sequence=fse_seq, name="fse-test")
    ft.tl.compute_ensemble(fd)
    ft.tl.compute_unpaired_probs(fd)

    tl_unpaired = fd.obs["unpaired_prob"].values

    # Compute using ViennaRNA directly
    fc = RNA.fold_compound(fse_seq)
    fc.pf()
    bpp = fc.bpp()
    n = len(fse_seq)

    vienna_unpaired = np.zeros(n)
    for i in range(1, n + 1):
        paired_prob = 0.0
        for j in range(1, n + 1):
            if i != j:
                if i < j:
                    paired_prob += bpp[i][j]
                else:
                    paired_prob += bpp[j][i]
        vienna_unpaired[i - 1] = 1.0 - paired_prob

    # Check they match
    assert np.allclose(tl_unpaired, vienna_unpaired, atol=1e-6), (
        f"ft.tl unpaired probabilities don't match ViennaRNA on FSE!\n"
        f"Max difference: {np.max(np.abs(tl_unpaired - vienna_unpaired))}"
    )


def test_unpaired_prob_tl_path_random():
    """Test ft.tl.compute_unpaired_probs path on random sequence matches ViennaRNA."""
    # Random 120-nt sequence
    np.random.seed(42)
    bases = ["A", "C", "G", "U"]
    seq = "".join(np.random.choice(bases, 120))

    # Create FoldData and compute ensemble
    fd = FoldData(sequence=seq, name="random-test")
    ft.tl.compute_ensemble(fd)
    ft.tl.compute_unpaired_probs(fd)

    tl_unpaired = fd.obs["unpaired_prob"].values

    # Compute using ViennaRNA directly
    fc = RNA.fold_compound(seq)
    fc.pf()
    bpp = fc.bpp()
    n = len(seq)

    vienna_unpaired = np.zeros(n)
    for i in range(1, n + 1):
        paired_prob = 0.0
        for j in range(1, n + 1):
            if i != j:
                if i < j:
                    paired_prob += bpp[i][j]
                else:
                    paired_prob += bpp[j][i]
        vienna_unpaired[i - 1] = 1.0 - paired_prob

    # Check they match
    assert np.allclose(tl_unpaired, vienna_unpaired, atol=1e-6), (
        f"ft.tl unpaired probabilities don't match ViennaRNA on random sequence!\n"
        f"Max difference: {np.max(np.abs(tl_unpaired - vienna_unpaired))}"
    )
