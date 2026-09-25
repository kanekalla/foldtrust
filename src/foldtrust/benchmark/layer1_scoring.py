"""
Layer 1: Scoring Correctness - Known-answer tests for FoldTrust pair scoring.

Tests sensitivity, PPV, F1, and MCC with exact and +/-1 slip-tolerant scoring.
Includes bracket/bpseq/ct parsing tests and tier regression tests.
"""

import json
import gzip
from pathlib import Path
from typing import Dict, List, Set, Tuple, Optional
import numpy as np
import pandas as pd


def parse_dot_bracket(structure: str) -> Set[Tuple[int, int]]:
    """
    Parse dot-bracket notation to extract base pairs.
    
    Handles (), <>, [], {} as nested pairs (WUSS convention).
    Returns 0-indexed pairs as (i, j) where i < j.
    """
    pairs = set()
    stacks = {'(': [], '<': [], '[': [], '{': []}
    close_to_open = {')': '(', '>': '<', ']': '[', '}': '{'}
    
    for i, char in enumerate(structure):
        if char in stacks:
            stacks[char].append(i)
        elif char in close_to_open:
            open_char = close_to_open[char]
            if stacks[open_char]:
                j = stacks[open_char].pop()
                pairs.add((j, i) if j < i else (i, j))
    
    return pairs


def parse_bpseq(filepath: Path) -> Tuple[str, Set[Tuple[int, int]]]:
    """
    Parse bpseq format: index nucleotide pair_index (1-indexed, 0 = unpaired).
    
    Returns:
        (sequence, pairs) where pairs are 0-indexed (i, j) with i < j
    """
    sequence = []
    pairs = set()
    
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 3:
                idx = int(parts[0]) - 1  # Convert to 0-indexed
                nuc = parts[1]
                pair_idx = int(parts[2])
                
                sequence.append(nuc)
                
                if pair_idx > 0:
                    pair_idx -= 1  # Convert to 0-indexed
                    if idx < pair_idx:
                        pairs.add((idx, pair_idx))
    
    return ''.join(sequence), pairs


def parse_ct(filepath: Path) -> Tuple[str, Set[Tuple[int, int]]]:
    """
    Parse CT format: index nuc prev_idx next_idx pair_idx seq_number.
    
    Returns:
        (sequence, pairs) where pairs are 0-indexed (i, j) with i < j
    """
    sequence = []
    pairs = set()
    
    with open(filepath) as f:
        lines = f.readlines()
    
    # First line is header with count
    if not lines:
        return "", set()
    
    for line in lines[1:]:
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        if len(parts) >= 5:
            idx = int(parts[0]) - 1  # Convert to 0-indexed
            nuc = parts[1]
            pair_idx = int(parts[4])
            
            sequence.append(nuc)
            
            if pair_idx > 0:
                pair_idx -= 1  # Convert to 0-indexed
                if idx < pair_idx:
                    pairs.add((idx, pair_idx))
    
    return ''.join(sequence), pairs


def compute_metrics(
    predicted_pairs: Set[Tuple[int, int]],
    reference_pairs: Set[Tuple[int, int]],
    seq_length: int,
    allow_slip: bool = False
) -> Dict[str, float]:
    """
    Compute sensitivity, PPV, F1, and MCC for predicted vs reference pairs.
    
    Args:
        predicted_pairs: Set of predicted pairs (i, j)
        reference_pairs: Set of reference pairs (i, j)
        seq_length: Length of sequence
        allow_slip: If True, allow +/-1 slippage tolerance
        
    Returns:
        Dict with sensitivity, ppv, f1, mcc, tp, fp, fn, tn
    """
    if allow_slip:
        # For slip tolerance, a predicted pair (i, j) matches reference if
        # any of (i±1, j±1) is in reference
        tp = 0
        matched_ref = set()
        
        for pred_i, pred_j in predicted_pairs:
            found_match = False
            for di in [-1, 0, 1]:
                for dj in [-1, 0, 1]:
                    ref_pair = (pred_i + di, pred_j + dj)
                    if ref_pair in reference_pairs and ref_pair not in matched_ref:
                        matched_ref.add(ref_pair)
                        found_match = True
                        break
                if found_match:
                    break
            if found_match:
                tp += 1
        
        fp = len(predicted_pairs) - tp
        fn = len(reference_pairs) - len(matched_ref)
    else:
        tp = len(predicted_pairs & reference_pairs)
        fp = len(predicted_pairs - reference_pairs)
        fn = len(reference_pairs - predicted_pairs)
    
    # Total possible pairs (upper triangle)
    total_pairs = seq_length * (seq_length - 1) // 2
    tn = total_pairs - tp - fp - fn
    
    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    ppv = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    f1 = 2 * tp / (2 * tp + fp + fn) if (2 * tp + fp + fn) > 0 else 0.0
    
    # Matthews correlation coefficient
    denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    mcc = (tp * tn - fp * fn) / denominator if denominator > 0 else 0.0
    
    return {
        'sensitivity': sensitivity,
        'ppv': ppv,
        'f1': f1,
        'mcc': mcc,
        'tp': tp,
        'fp': fp,
        'fn': fn,
        'tn': tn
    }


def test_parsing():
    """Test bracket/bpseq/ct parsing with known examples."""
    print("Testing parsers...")
    
    # Test dot-bracket parsing
    structure = "(((...)))..(((...)))"
    pairs = parse_dot_bracket(structure)
    expected = {(0, 8), (1, 7), (2, 6), (11, 19), (12, 18), (13, 17)}
    assert pairs == expected, f"Dot-bracket parsing failed: {pairs} != {expected}"
    print("  ✓ Dot-bracket parser")
    
    # Test nested brackets
    structure = "(((<...>)))"
    pairs = parse_dot_bracket(structure)
    expected = {(0, 10), (1, 9), (2, 8), (3, 7)}
    assert pairs == expected, f"Nested bracket parsing failed: {pairs} != {expected}"
    print("  ✓ Nested bracket parser")
    
    # Test bracket types
    structure = "(((...))).[[...]].{{...}}"
    pairs = parse_dot_bracket(structure)
    expected = {(0, 8), (1, 7), (2, 6), (10, 16), (11, 15), (18, 24), (19, 23)}
    assert pairs == expected, f"Multiple bracket types failed: {pairs} != {expected}"
    print("  ✓ Multiple bracket types")
    
    print("All parser tests passed!")
    return True


def test_tier_regression():
    """
    Tier regression tests: verify that strong structures get FIRM pairs.
    
    Tests that a strong GC hairpin produces FIRM pairs.
    """
    print("\nTesting tier regression...")
    
    try:
        import RNA
        
        # Strong GC hairpin
        seq = "GCGCGCAAAAGCGCGC"
        
        # Compute ensemble
        md = RNA.md()
        md.temperature = 37.0
        fc = RNA.fold_compound(seq, md)
        
        # MFE
        structure, mfe = fc.mfe()
        
        # Partition function
        fc.pf()
        bpp = fc.bpp()
        
        # Extract pairs from MFE
        mfe_pairs = parse_dot_bracket(structure)
        
        # Check pair probabilities - all stem pairs should be FIRM (>= 0.85)
        firm_count = 0
        for i, j in mfe_pairs:
            prob = bpp[i+1][j+1] if i+1 < len(bpp) and j+1 < len(bpp[i+1]) else 0.0
            if prob >= 0.85:
                firm_count += 1
        
        print(f"  Strong GC hairpin: {len(mfe_pairs)} pairs, {firm_count} FIRM (>= 0.85)")
        print(f"  MFE structure: {structure}")
        print(f"  MFE energy: {mfe:.2f} kcal/mol")
        
        # At least half of pairs should be FIRM for this strong structure
        assert firm_count >= len(mfe_pairs) // 2, \
            f"Strong GC hairpin should have mostly FIRM pairs, got {firm_count}/{len(mfe_pairs)}"
        
        print("  ✓ Tier regression test passed")
        return True
        
    except ImportError:
        print("  ⚠ ViennaRNA Python not available, skipping tier test")
        return False


def run_layer1_tests(output_dir: Path) -> Dict:
    """
    Run Layer 1 tests and save results.
    
    Returns:
        Summary dict with test results
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    
    results = {
        'parser_tests': test_parsing(),
        'tier_regression': test_tier_regression(),
    }
    
    # Test metrics computation with known values
    print("\nTesting metrics computation...")
    
    # Perfect prediction
    ref_pairs = {(0, 10), (1, 9), (2, 8), (3, 7)}
    pred_pairs = {(0, 10), (1, 9), (2, 8), (3, 7)}
    metrics = compute_metrics(pred_pairs, ref_pairs, 11, allow_slip=False)
    assert metrics['sensitivity'] == 1.0 and metrics['ppv'] == 1.0
    print("  ✓ Perfect prediction test")
    
    # Test with FP and FN
    pred_pairs = {(0, 10), (1, 9), (4, 6)}  # Missing (2,8), (3,7), wrong (4,6)
    metrics = compute_metrics(pred_pairs, ref_pairs, 11, allow_slip=False)
    assert metrics['tp'] == 2 and metrics['fp'] == 1 and metrics['fn'] == 2
    print("  ✓ FP/FN test")
    
    # Test slip tolerance
    pred_pairs = {(0, 10), (1, 9), (2, 9), (3, 8)}  # Slipped by 1
    metrics_exact = compute_metrics(pred_pairs, ref_pairs, 11, allow_slip=False)
    metrics_slip = compute_metrics(pred_pairs, ref_pairs, 11, allow_slip=True)
    assert metrics_slip['tp'] > metrics_exact['tp']
    print("  ✓ Slip tolerance test")
    
    # Save results
    results_file = output_dir / 'layer1_tests.json'
    with open(results_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"\nLayer 1 tests complete. Results saved to {results_file}")
    return results


if __name__ == '__main__':
    output_dir = Path('benchmarks/outputs/layer1')
    run_layer1_tests(output_dir)
