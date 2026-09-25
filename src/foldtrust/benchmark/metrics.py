"""Core metrics for structure prediction and calibration evaluation."""

from typing import Dict, List, Set, Tuple

import numpy as np
from sklearn.metrics import average_precision_score, roc_auc_score


def parse_structure_pairs(structure: str) -> set:
    """Parse base pairs from dot-bracket notation."""
    pairs = set()
    stack = []

    for i, char in enumerate(structure):
        if char in "([{<":
            stack.append(i)
        elif char in ")]}>" and stack:
            j = stack.pop()
            pairs.add((j, i) if j < i else (i, j))

    return pairs


def compute_accuracy_metrics(
    predicted: str, reference: str, allow_slip: bool = False
) -> Dict[str, float]:
    """
    Compute sensitivity, PPV, F1, and MCC for predicted vs reference structure.

    Args:
        predicted: Predicted structure in dot-bracket notation
        reference: Reference structure in dot-bracket notation
        allow_slip: If True, allow one-nucleotide slippage

    Returns:
        Dictionary with sensitivity, ppv, f1, and mcc
    """
    pred_pairs = parse_structure_pairs(predicted)
    ref_pairs = parse_structure_pairs(reference)

    if allow_slip:
        tp = 0
        for pred_pair in pred_pairs:
            i, j = pred_pair
            if pred_pair in ref_pairs:
                tp += 1
            elif any(
                (i + di, j + dj) in ref_pairs
                for di in [-1, 0, 1]
                for dj in [-1, 0, 1]
                if not (di == 0 and dj == 0)
            ):
                tp += 1
    else:
        tp = len(pred_pairs & ref_pairs)

    fp = len(pred_pairs) - tp
    fn = len(ref_pairs) - tp
    tn = len(reference) * (len(reference) - 1) // 2 - tp - fp - fn

    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    ppv = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    f1 = 2 * tp / (2 * tp + fp + fn) if (2 * tp + fp + fn) > 0 else 0.0

    denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    mcc = (tp * tn - fp * fn) / denominator if denominator > 0 else 0.0

    return {
        "sensitivity": sensitivity,
        "ppv": ppv,
        "f1": f1,
        "mcc": mcc,
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "tn": tn,
    }


def compute_calibration_metrics(
    pair_probs: np.ndarray, reference_pairs: set, n_bins: int = 10
) -> Dict:
    """
    Compute calibration metrics: reliability diagram data, ECE, AUROC, AUPRC.

    Args:
        pair_probs: NxN matrix of predicted pair probabilities
        reference_pairs: Set of true base pairs (i, j) tuples
        n_bins: Number of bins for calibration curve

    Returns:
        Dictionary with bin_means, bin_accs, ece, auroc, auprc
    """
    n = pair_probs.shape[0]

    y_true = []
    y_prob = []

    for i in range(n):
        for j in range(i + 1, n):
            y_true.append(1 if (i, j) in reference_pairs else 0)
            y_prob.append(pair_probs[i, j])

    y_true = np.array(y_true)
    y_prob = np.array(y_prob)

    if len(y_true) == 0 or np.sum(y_true) == 0:
        return {
            "bin_means": [],
            "bin_accs": [],
            "bin_counts": [],
            "ece": 0.0,
            "auroc": 0.0,
            "auprc": 0.0,
        }

    auroc = roc_auc_score(y_true, y_prob) if np.sum(y_true) > 0 else 0.0
    auprc = average_precision_score(y_true, y_prob) if np.sum(y_true) > 0 else 0.0

    bin_edges = np.linspace(0, 1, n_bins + 1)
    bin_means = []
    bin_accs = []
    bin_counts = []
    ece = 0.0

    for b in range(n_bins):
        mask = (y_prob >= bin_edges[b]) & (y_prob < bin_edges[b + 1])
        if b == n_bins - 1:
            mask = mask | (y_prob == 1.0)

        if np.sum(mask) > 0:
            bin_mean = np.mean(y_prob[mask])
            bin_acc = np.mean(y_true[mask])
            bin_count = np.sum(mask)

            bin_means.append(bin_mean)
            bin_accs.append(bin_acc)
            bin_counts.append(int(bin_count))

            ece += (bin_count / len(y_true)) * abs(bin_mean - bin_acc)

    return {
        "bin_means": bin_means,
        "bin_accs": bin_accs,
        "bin_counts": bin_counts,
        "ece": ece,
        "auroc": auroc,
        "auprc": auprc,
    }


def compute_probing_metrics(
    unpaired_probs: np.ndarray, reactivities: np.ndarray
) -> Dict[str, float]:
    """
    Compute correlation between unpaired probabilities and experimental reactivity.

    Args:
        unpaired_probs: Array of unpaired probabilities per nucleotide
        reactivities: Array of experimental reactivity values (e.g., SHAPE)

    Returns:
        Dictionary with spearman correlation and AUROC if binary classification possible
    """
    from scipy.stats import spearmanr

    valid_mask = ~np.isnan(reactivities)
    if np.sum(valid_mask) == 0:
        return {"spearman": 0.0, "auroc": 0.0}

    unpaired_valid = unpaired_probs[valid_mask]
    react_valid = reactivities[valid_mask]

    spearman, p_value = spearmanr(unpaired_valid, react_valid)

    thresh = np.median(react_valid)
    binary_labels = (react_valid > thresh).astype(int)

    auroc = 0.0
    if len(np.unique(binary_labels)) > 1:
        auroc = roc_auc_score(binary_labels, unpaired_valid)

    return {
        "spearman": spearman,
        "spearman_pvalue": p_value,
        "auroc": auroc,
    }


def stratify_by_length(sequences: List[str], bins: List[int]) -> Dict[str, List[int]]:
    """
    Stratify sequence indices by length bins.

    Args:
        sequences: List of RNA sequences
        bins: List of bin edges (e.g., [0, 50, 100, 200, 500])

    Returns:
        Dictionary mapping bin labels to lists of indices
    """
    stratified = {}

    for i, seq in enumerate(sequences):
        length = len(seq)

        for b in range(len(bins) - 1):
            if bins[b] <= length < bins[b + 1]:
                label = f"{bins[b]}-{bins[b + 1]}"
                if label not in stratified:
                    stratified[label] = []
                stratified[label].append(i)
                break
        else:
            if length >= bins[-1]:
                label = f"{bins[-1]}+"
                if label not in stratified:
                    stratified[label] = []
                stratified[label].append(i)

    return stratified


def compute_tier_accuracy(stems: List[Dict], reference_pairs: set) -> Dict[str, Dict[str, float]]:
    """
    Compute PPV for pairs in each reliability tier (firm/soft/floppy).

    Args:
        stems: List of stem dictionaries with 'flag' and 'pairs' keys
        reference_pairs: Set of true base pairs

    Returns:
        Dictionary mapping tier names to accuracy metrics
    """
    tier_stats = {
        "firm": {"tp": 0, "total": 0},
        "soft": {"tp": 0, "total": 0},
        "floppy": {"tp": 0, "total": 0},
    }

    for stem in stems:
        tier = stem["flag"]
        for pair in stem["pairs"]:
            tier_stats[tier]["total"] += 1
            if pair in reference_pairs:
                tier_stats[tier]["tp"] += 1

    tier_accuracy = {}
    for tier, stats in tier_stats.items():
        if stats["total"] > 0:
            tier_accuracy[tier] = {
                "ppv": stats["tp"] / stats["total"],
                "count": stats["total"],
            }
        else:
            tier_accuracy[tier] = {"ppv": 0.0, "count": 0}

    return tier_accuracy


def remove_pseudoknots(pairs: Set[Tuple[int, int]]) -> Set[Tuple[int, int]]:
    """
    Remove pseudoknots using greedy algorithm (keep longest compatible pairs).

    Args:
        pairs: Set of base pairs (i, j) tuples

    Returns:
        Set of non-crossing base pairs
    """
    if not pairs:
        return set()

    # Sort by span length (descending)
    sorted_pairs = sorted(pairs, key=lambda p: p[1] - p[0], reverse=True)

    canonical = []
    for pair in sorted_pairs:
        # Check if this pair crosses any accepted pair
        crosses = False
        for accepted in canonical:
            a1, a2 = accepted
            p1, p2 = pair
            # Crossing condition: a1 < p1 < a2 < p2 or p1 < a1 < p2 < a2
            if (a1 < p1 < a2 < p2) or (p1 < a1 < p2 < a2):
                crosses = True
                break

        if not crosses:
            canonical.append(pair)

    return set(canonical)


def compute_structure_metrics(
    ref_pairs: Set[Tuple[int, int]], pred_pairs: Set[Tuple[int, int]]
) -> Dict[str, float]:
    """
    Compute exact structure accuracy metrics.

    Args:
        ref_pairs: Set of reference base pairs
        pred_pairs: Set of predicted base pairs

    Returns:
        Dictionary with sensitivity, ppv, f1, tp, fp, fn
    """
    tp = len(ref_pairs & pred_pairs)
    fp = len(pred_pairs - ref_pairs)
    fn = len(ref_pairs - pred_pairs)

    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    ppv = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    f1 = 2 * tp / (2 * tp + fp + fn) if (2 * tp + fp + fn) > 0 else 0.0

    return {
        "sensitivity": sensitivity,
        "ppv": ppv,
        "f1": f1,
        "tp": tp,
        "fp": fp,
        "fn": fn,
    }


def compute_slip_tolerant_metrics(
    ref_pairs: Set[Tuple[int, int]], pred_pairs: Set[Tuple[int, int]]
) -> Dict[str, float]:
    """
    Compute slip-tolerant structure accuracy metrics.

    A predicted pair (i,j) is a slip-TP if any of (i,j), (i±1,j), (i,j±1) is in reference.
    A reference pair is recovered if any of (i,j), (i±1,j), (i,j±1) is predicted.

    Args:
        ref_pairs: Set of reference base pairs
        pred_pairs: Set of predicted base pairs

    Returns:
        Dictionary with sensitivity, ppv, f1
    """
    # Count slip-TP for PPV: predicted pairs with a match
    slip_tp_pred = 0
    for i, j in pred_pairs:
        neighbors = {(i, j), (i - 1, j), (i + 1, j), (i, j - 1), (i, j + 1)}
        if neighbors & ref_pairs:
            slip_tp_pred += 1

    # Count recovered reference pairs for sensitivity
    slip_tp_ref = 0
    for i, j in ref_pairs:
        neighbors = {(i, j), (i - 1, j), (i + 1, j), (i, j - 1), (i, j + 1)}
        if neighbors & pred_pairs:
            slip_tp_ref += 1

    ppv = slip_tp_pred / len(pred_pairs) if len(pred_pairs) > 0 else 0.0
    sensitivity = slip_tp_ref / len(ref_pairs) if len(ref_pairs) > 0 else 0.0
    f1 = 2 * ppv * sensitivity / (ppv + sensitivity) if (ppv + sensitivity) > 0 else 0.0

    return {
        "sensitivity": sensitivity,
        "ppv": ppv,
        "f1": f1,
    }
