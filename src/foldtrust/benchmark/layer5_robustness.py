"""Layer 5: Robustness analysis across temperature, parameter sets, and window context.

This module tests FoldTrust's structure and tier predictions under three types of perturbations:
1. Temperature variations (25, 30, 37, 42°C; 37°C is the reference)
2. Energy parameter sets (Turner2004 reference vs Andronescu2007 vs Langdon2018)
3. Window context: extend by real flanking sequence (0, 25, 50, 100 nt on each side)

All conditions are compared to the 37°C Turner2004 baseline for the core window.
Retention metrics are computed only over stems/pairs that exist in the reference condition.
"""

import json
import urllib.error
import urllib.request
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

try:
    import RNA

    HAS_RNA = True
except ImportError:
    HAS_RNA = False

from foldtrust.utils import read_fasta
from foldtrust.vienna import (
    ViennaRNAError,
    compute_pair_probs_with_params,
    fold_with_params,
    parse_stems,
)

# Verified case coordinates from Layer 0
CASE_COORDS = {
    "sars2-fse": {
        "accession": "NC_045512.2",
        "start": 13462,
        "end": 13542,
        "length": 81,
    },
    "smn2-iss-n1": {
        "accession": "NG_008728.1",
        "start": 31999,
        "end": 32152,
        "length": 154,
    },
    "cftr-5utr": {
        "accession": "NM_000492.4",
        "start": 1,
        "end": 200,
        "length": 200,
    },
    "mapt-e10": {
        "accession": "NG_007398.2",
        "start": 120818,
        "end": 121000,
        "length": 183,
    },
    "hcv-ires-dii": {
        "accession": "AF009606.1",
        "start": 44,
        "end": 118,
        "length": 75,
    },
}


def fetch_fasta_from_ncbi(accession: str, cache_dir: Path) -> Optional[str]:
    """
    Fetch complete FASTA record from NCBI, with caching.

    Args:
        accession: GenBank/RefSeq accession
        cache_dir: Directory for cached files

    Returns:
        FASTA content string, or None if fetch fails
    """
    cache_file = cache_dir / f"{accession}.fa"

    if cache_file.exists():
        return cache_file.read_text()

    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={accession}&rettype=fasta&retmode=text"

    try:
        with urllib.request.urlopen(url, timeout=30) as response:
            content = response.read().decode("utf-8")

        cache_dir.mkdir(parents=True, exist_ok=True)
        cache_file.write_text(content)

        return content
    except urllib.error.URLError as e:
        print(f"  Warning: Failed to fetch {accession}: {e}")
        return None


def extract_sequence_from_fasta(fasta_content: str) -> str:
    """Extract sequence from FASTA content (strip header)."""
    lines = fasta_content.strip().split("\n")
    return "".join(line.strip() for line in lines if not line.startswith(">"))


def get_flanking_sequences(
    case_name: str,
    flank_size: int,
    cache_dir: Path,
) -> Tuple[str, str]:
    """
    Fetch real flanking sequences from NCBI.

    Args:
        case_name: Case identifier
        flank_size: Number of nucleotides to fetch on each side
        cache_dir: Cache directory

    Returns:
        (left_flank, right_flank) tuple of sequences
    """
    coords = CASE_COORDS[case_name]
    accession = coords["accession"]
    start = coords["start"]
    end = coords["end"]

    # Fetch full record
    fasta_content = fetch_fasta_from_ncbi(accession, cache_dir)
    if not fasta_content:
        return "", ""

    full_seq = extract_sequence_from_fasta(fasta_content).upper().replace("T", "U")

    # Extract left flank (clip at record start)
    left_start = max(0, start - 1 - flank_size)
    left_end = start - 1
    left_flank = full_seq[left_start:left_end]

    # Extract right flank (clip at record end)
    right_start = end
    right_end = min(len(full_seq), end + flank_size)
    right_flank = full_seq[right_start:right_end]

    return left_flank, right_flank


def verify_core_sequence(
    case_name: str,
    core_sequence: str,
    cache_dir: Path,
) -> bool:
    """
    Verify that core sequence is an exact substring of the cached record.

    Returns:
        True if verified
    """
    coords = CASE_COORDS[case_name]
    accession = coords["accession"]
    start = coords["start"]
    end = coords["end"]

    fasta_content = fetch_fasta_from_ncbi(accession, cache_dir)
    if not fasta_content:
        return False

    full_seq = extract_sequence_from_fasta(fasta_content).upper().replace("T", "U")
    expected_seq = full_seq[start - 1 : end]

    core_norm = core_sequence.upper().replace("T", "U")

    return core_norm == expected_seq


def compute_mea_structure(
    sequence: str,
    param_set: str = "Turner2004",
    temperature: float = 37.0,
    gamma: float = 1.0,
) -> Tuple[List[Tuple[int, int]], float]:
    """
    Compute MEA structure using ViennaRNA's fc.MEA(gamma).

    Args:
        sequence: RNA sequence
        param_set: Parameter set name
        temperature: Temperature in Celsius
        gamma: MEA gamma parameter (default 1.0)

    Returns:
        (pairs, mea_energy) where pairs is list of (i, j) 0-based tuples
    """
    if not HAS_RNA:
        raise ViennaRNAError("ViennaRNA not available")

    # Use subprocess like fold_with_params to ensure parameter isolation
    import subprocess
    import sys

    code = f"""
import RNA

sequence = {sequence!r}
param_set = {param_set!r}
temperature = {temperature}
gamma = {gamma}

# Load parameters
if param_set == "Turner2004":
    RNA.params_load_RNA_Turner2004()
elif param_set == "Andronescu2007":
    RNA.params_load_RNA_Andronescu2007()
elif param_set == "Langdon2018":
    RNA.params_load_RNA_Langdon2018()

# Create fold compound
md = RNA.md()
md.temperature = temperature
fc = RNA.fold_compound(sequence, md)

# Compute partition function
fc.pf()

# Compute MEA
mea_result = fc.MEA(gamma)
mea_struct = mea_result[0]
mea_energy = mea_result[1]

print(f"{{mea_struct}}|{{mea_energy}}")
"""

    try:
        result = subprocess.run(
            [sys.executable, "-c", code],
            capture_output=True,
            text=True,
            timeout=30,
            check=True,
        )

        output = result.stdout.strip()
        mea_struct, mea_energy_str = output.split("|")
        mea_energy = float(mea_energy_str)

        # Parse structure to pairs
        pairs = []
        stack = []
        for i, char in enumerate(mea_struct):
            if char == "(":
                stack.append(i)
            elif char == ")" and stack:
                j = stack.pop()
                pairs.append((j, i))

        return pairs, mea_energy

    except subprocess.CalledProcessError as e:
        raise ViennaRNAError(f"MEA computation failed: {e.stderr}")
    except Exception as e:
        raise ViennaRNAError(f"MEA error: {e}")


def compute_ensemble_defect(
    sequence: str,
    structure: str,
    param_set: str = "Turner2004",
    temperature: float = 37.0,
) -> float:
    """
    Compute ensemble defect: expected number of incorrectly paired bases.

    Ensemble defect = sum over positions of (1 - probability of being in correct state)

    Args:
        sequence: RNA sequence
        structure: Dot-bracket structure
        param_set: Parameter set
        temperature: Temperature

    Returns:
        Ensemble defect
    """
    prob_matrix = compute_pair_probs_with_params(sequence, param_set, temperature)
    n = len(sequence)
    defect = 0.0

    # Parse structure to get paired positions
    paired = {}
    stack = []
    for i, char in enumerate(structure):
        if char == "(":
            stack.append(i)
        elif char == ")" and stack:
            j = stack.pop()
            paired[j] = i
            paired[i] = j

    # For each position, probability of correct state
    for i in range(n):
        if i in paired:
            j = paired[i]
            correct_prob = prob_matrix[i, j]
        else:
            correct_prob = 1.0 - np.sum(prob_matrix[i, :])

        defect += 1.0 - correct_prob

    return defect


def compute_basepair_distance(pairs1: List[Tuple[int, int]], pairs2: List[Tuple[int, int]]) -> int:
    """
    Compute base-pair distance: number of pairs in symmetric difference.

    Args:
        pairs1, pairs2: Lists of (i, j) pairs

    Returns:
        |pairs1 - pairs2| + |pairs2 - pairs1|
    """
    set1 = set(pairs1)
    set2 = set(pairs2)
    return len(set1 - set2) + len(set2 - set1)


def compute_stem_retention(
    ref_stems: List[Dict],
    cond_prob_matrix: np.ndarray,
    cond_mea_pairs: List[Tuple[int, int]],
) -> Dict:
    """
    Compute retention of reference stems in a condition.

    For each reference stem tier (FIRM/SOFT/FLOPPY), compute:
    - Fraction of stems whose pairs are still paired in condition MEA
    - Mean change in stem mean probability

    Args:
        ref_stems: Reference stems from parse_stems
        cond_prob_matrix: Condition probability matrix
        cond_mea_pairs: Condition MEA pairs

    Returns:
        Dict with retention metrics per tier
    """
    cond_pairs_set = set(cond_mea_pairs)

    # Group stems by tier
    tiers = {"firm": [], "soft": [], "floppy": []}
    for stem in ref_stems:
        tiers[stem["flag"]].append(stem)

    results = {}
    for tier_name, stems in tiers.items():
        if not stems:
            results[tier_name] = {
                "count": 0,
                "retention": None,
                "mean_prob_change": None,
            }
            continue

        # Count how many stems are retained (all pairs still paired in MEA)
        retained = 0
        prob_changes = []

        for stem in stems:
            stem_pairs = stem["pairs"]
            ref_mean_prob = stem["mean_prob"]

            # Check if all pairs in stem are paired in condition MEA
            all_paired = all(pair in cond_pairs_set for pair in stem_pairs)
            if all_paired:
                retained += 1

            # Compute mean probability change
            cond_probs = [cond_prob_matrix[i, j] for i, j in stem_pairs]
            cond_mean_prob = np.mean(cond_probs)
            prob_changes.append(cond_mean_prob - ref_mean_prob)

        results[tier_name] = {
            "count": len(stems),
            "retention": retained / len(stems),
            "mean_prob_change": np.mean(prob_changes),
        }

    return results


def compute_baseline_retention(
    ref_stems: List[Dict],
    ref_mea_pairs: List[Tuple[int, int]],
) -> Dict:
    """
    Compute baseline retention: MFE stems retained in MEA.

    This measures the MFE-vs-MEA gap for the same condition (37°C Turner2004).

    Args:
        ref_stems: Stems from MFE structure
        ref_mea_pairs: MEA pairs from same condition

    Returns:
        Retention dict by tier
    """
    mea_pairs_set = set(ref_mea_pairs)

    # Group stems by tier
    tiers = {"firm": [], "soft": [], "floppy": []}
    for stem in ref_stems:
        tiers[stem["flag"]].append(stem)

    results = {}
    for tier_name, stems in tiers.items():
        if not stems:
            results[tier_name] = {
                "count": 0,
                "retention": None,
            }
            continue

        # Count stems whose pairs are all in MEA
        retained = 0
        for stem in stems:
            all_paired = all(pair in mea_pairs_set for pair in stem["pairs"])
            if all_paired:
                retained += 1

        results[tier_name] = {
            "count": len(stems),
            "retention": retained / len(stems),
        }

    return results


def analyze_temperature_condition(
    sequence: str,
    temp: float,
    ref_structure: str,
    ref_mea_pairs: List[Tuple[int, int]],
    ref_stems: List[Dict],
) -> Dict:
    """
    Analyze one temperature condition against 37°C Turner2004 reference.

    Returns dict with metrics: mfe_energy, ensemble_defect, bp_distance_mfe,
    bp_distance_mea, stem_retention (by tier)
    """
    param_set = "Turner2004"

    # Compute MFE
    structure, mfe_energy = fold_with_params(sequence, param_set, temp)

    # Parse MFE pairs
    mfe_pairs = []
    stack = []
    for i, char in enumerate(structure):
        if char == "(":
            stack.append(i)
        elif char == ")" and stack:
            j = stack.pop()
            mfe_pairs.append((j, i))

    # Compute MEA
    mea_pairs, mea_energy = compute_mea_structure(sequence, param_set, temp)

    # Ensemble defect
    ensemble_defect = compute_ensemble_defect(sequence, structure, param_set, temp)

    # Base-pair distances
    bp_dist_mfe = compute_basepair_distance(
        [(i, j) for i, j in mfe_pairs], [(i, j) for i, j in ref_mea_pairs]
    )
    bp_dist_mea = compute_basepair_distance(mea_pairs, ref_mea_pairs)

    # Probability matrix for stem retention
    prob_matrix = compute_pair_probs_with_params(sequence, param_set, temp)

    # Stem retention
    stem_retention = compute_stem_retention(ref_stems, prob_matrix, mea_pairs)

    return {
        "mfe_energy": mfe_energy,
        "ensemble_defect": ensemble_defect,
        "bp_distance_mfe": bp_dist_mfe,
        "bp_distance_mea": bp_dist_mea,
        "stem_retention": stem_retention,
    }


def analyze_params_condition(
    sequence: str,
    param_set: str,
    ref_structure: str,
    ref_mea_pairs: List[Tuple[int, int]],
    ref_stems: List[Dict],
) -> Dict:
    """
    Analyze one parameter set against Turner2004 reference.
    """
    temp = 37.0

    # Compute MFE
    structure, mfe_energy = fold_with_params(sequence, param_set, temp)

    # Parse MFE pairs
    mfe_pairs = []
    stack = []
    for i, char in enumerate(structure):
        if char == "(":
            stack.append(i)
        elif char == ")" and stack:
            j = stack.pop()
            mfe_pairs.append((j, i))

    # Compute MEA
    mea_pairs, mea_energy = compute_mea_structure(sequence, param_set, temp)

    # Ensemble defect
    ensemble_defect = compute_ensemble_defect(sequence, structure, param_set, temp)

    # Base-pair distances
    bp_dist_mfe = compute_basepair_distance(mfe_pairs, ref_mea_pairs)
    bp_dist_mea = compute_basepair_distance(mea_pairs, ref_mea_pairs)

    # Probability matrix
    prob_matrix = compute_pair_probs_with_params(sequence, param_set, temp)

    # Stem retention
    stem_retention = compute_stem_retention(ref_stems, prob_matrix, mea_pairs)

    return {
        "mfe_energy": mfe_energy,
        "ensemble_defect": ensemble_defect,
        "bp_distance_mfe": bp_dist_mfe,
        "bp_distance_mea": bp_dist_mea,
        "stem_retention": stem_retention,
    }


def analyze_window_context(
    core_sequence: str,
    left_flank: str,
    right_flank: str,
    flank_size: int,
    case_name: str,
    ref_structure: str,
    ref_mea_pairs: List[Tuple[int, int]],
    ref_stems: List[Dict],
) -> Dict:
    """
    Analyze window with flanking context.

    Tests retention of core window stems when flanking sequence is added.
    Only the core window stems (from ref_stems) are counted for retention.
    BP distance counts only pairs with both ends in the core window.

    Args:
        core_sequence: Core window sequence
        left_flank: Left flanking sequence
        right_flank: Right flanking sequence
        flank_size: Requested flank size (for metadata)
        case_name: Case name for coordinate metadata
        ref_structure: Reference structure (for core window)
        ref_mea_pairs: Reference MEA pairs (0-based, core window coords)
        ref_stems: Reference stems (from parse_stems on core window)

    Returns:
        Dict with metrics
    """
    # Construct extended sequence
    extended_seq = left_flank + core_sequence + right_flank
    core_offset = len(left_flank)
    core_end = core_offset + len(core_sequence)

    # Fold extended sequence
    param_set = "Turner2004"
    temp = 37.0

    structure, mfe_energy = fold_with_params(extended_seq, param_set, temp)
    mea_pairs, mea_energy = compute_mea_structure(extended_seq, param_set, temp)

    # Map core window reference pairs to extended coordinates
    ref_pairs_extended = [(i + core_offset, j + core_offset) for i, j in ref_mea_pairs]

    # Filter MEA pairs to only those with BOTH ends in core window
    mea_pairs_core_only = [
        (i, j) for i, j in mea_pairs if core_offset <= i < core_end and core_offset <= j < core_end
    ]

    # Base-pair distance (core window pairs only)
    bp_dist_mea = compute_basepair_distance(ref_pairs_extended, mea_pairs_core_only)

    # Stem retention: check if core stems are retained in extended MEA
    # Need to translate stem pairs to extended coordinates
    extended_stems = []
    for stem in ref_stems:
        extended_pairs = [(i + core_offset, j + core_offset) for i, j in stem["pairs"]]
        extended_stem = {
            "id": stem["id"],
            "pairs": extended_pairs,
            "length": stem["length"],
            "mean_prob": stem["mean_prob"],
            "flag": stem["flag"],
        }
        extended_stems.append(extended_stem)

    # Compute retention using extended stems
    mea_pairs_set = set(mea_pairs)
    tier_retention = {"firm": [], "soft": [], "floppy": []}

    for stem in extended_stems:
        tier = stem["flag"]
        all_paired = all(pair in mea_pairs_set for pair in stem["pairs"])
        tier_retention[tier].append(1 if all_paired else 0)

    # Average retention per tier
    retention_results = {}
    for tier in ["firm", "soft", "floppy"]:
        if tier_retention[tier]:
            retention_results[tier] = {
                "count": len(tier_retention[tier]),
                "retention": sum(tier_retention[tier]) / len(tier_retention[tier]),
            }
        else:
            retention_results[tier] = {
                "count": 0,
                "retention": None,
            }

    # Compute flank coordinates and clipping status
    coords = CASE_COORDS[case_name]
    accession = coords["accession"]
    core_start = coords["start"]
    core_end = coords["end"]

    left_start = core_start - len(left_flank)
    left_end = core_start - 1
    left_clipped = len(left_flank) < flank_size

    right_start = core_end + 1
    right_end = core_end + len(right_flank)
    right_clipped = len(right_flank) < flank_size

    left_coords = f"{accession}:{left_start}-{left_end}" if len(left_flank) > 0 else None
    right_coords = f"{accession}:{right_start}-{right_end}" if len(right_flank) > 0 else None

    return {
        "flank_size": flank_size,
        "left_flank_len": len(left_flank),
        "right_flank_len": len(right_flank),
        "extended_len": len(extended_seq),
        "mfe_energy": mfe_energy,
        "bp_distance_mea": bp_dist_mea,
        "stem_retention": retention_results,
        "left_accession_coords": left_coords,
        "right_accession_coords": right_coords,
        "left_clipped": left_clipped,
        "right_clipped": right_clipped,
    }


def run_layer5_analysis(
    cases_dir: Path,
    output_dir: Path,
    cache_dir: Optional[Path] = None,
) -> Dict:
    """
    Run Layer 5 robustness analysis for all cases.

    For each case:
    1. Compute 37°C Turner2004 baseline
    2. Test temperatures: 25, 30, 42°C
    3. Test parameter sets: Andronescu2007, Langdon2018
    4. Test window context: extend by 0, 25, 50, 100 nt flanks

    Args:
        cases_dir: Path to data/cases
        output_dir: Path to benchmarks/outputs/layer5
        cache_dir: Cache directory for NCBI fetches (default: data/_cache)

    Returns:
        Summary dict
    """
    if not HAS_RNA:
        raise RuntimeError("ViennaRNA not available. Install: pip install ViennaRNA")

    if cache_dir is None:
        cache_dir = Path("data/_cache")

    output_dir.mkdir(parents=True, exist_ok=True)
    cache_dir.mkdir(parents=True, exist_ok=True)

    results = {}

    for case_name in CASE_COORDS.keys():
        print(f"\n=== {case_name} ===")

        case_dir = cases_dir / case_name
        seq_file = case_dir / "sequence.fa"

        if not seq_file.exists():
            print("  Skipping: no sequence.fa")
            continue

        _, sequence = read_fasta(seq_file)
        print(f"  Length: {len(sequence)} nt")

        # Verify sequence against NCBI
        print("  Verifying sequence...")
        verified = verify_core_sequence(case_name, sequence, cache_dir)
        if not verified:
            raise RuntimeError(
                f"{case_name}: Core sequence verification failed. "
                f"Sequence is not an exact substring of {CASE_COORDS[case_name]['accession']} "
                f"at coordinates {CASE_COORDS[case_name]['start']}-{CASE_COORDS[case_name]['end']}"
            )
        print("    ✓ Verified as exact substring of cached record")

        # Baseline: 37°C Turner2004
        print("  Computing baseline (37°C Turner2004)...")
        ref_structure, ref_mfe_energy = fold_with_params(sequence, "Turner2004", 37.0)
        ref_prob_matrix = compute_pair_probs_with_params(sequence, "Turner2004", 37.0)
        ref_stems = parse_stems(ref_structure, ref_prob_matrix)
        ref_mea_pairs, ref_mea_energy = compute_mea_structure(sequence, "Turner2004", 37.0)
        ref_ensemble_defect = compute_ensemble_defect(sequence, ref_structure, "Turner2004", 37.0)

        # Compute baseline retention (MFE stems in MEA for same condition)
        baseline_retention = compute_baseline_retention(ref_stems, ref_mea_pairs)

        case_results = {
            "case": case_name,
            "length": len(sequence),
            "accession": CASE_COORDS[case_name]["accession"],
            "coords": f"{CASE_COORDS[case_name]['start']}-{CASE_COORDS[case_name]['end']}",
            "baseline": {
                "mfe_energy": ref_mfe_energy,
                "ensemble_defect": ref_ensemble_defect,
                "n_stems": len(ref_stems),
                "stem_counts": {
                    "firm": sum(1 for s in ref_stems if s["flag"] == "firm"),
                    "soft": sum(1 for s in ref_stems if s["flag"] == "soft"),
                    "floppy": sum(1 for s in ref_stems if s["flag"] == "floppy"),
                },
                "retention": baseline_retention,  # MFE-in-MEA retention
            },
            "temperature": {},
            "parameters": {},
            "window_context": {},
        }

        # Add baseline row to temperature sweep (37°C)
        case_results["temperature"]["37C_baseline"] = {
            "mfe_energy": ref_mfe_energy,
            "ensemble_defect": ref_ensemble_defect,
            "bp_distance_mfe": 0,
            "bp_distance_mea": 0,
            "stem_retention": baseline_retention,
        }

        # Temperature sweep
        for temp in [25.0, 30.0, 42.0]:
            print(f"  Temperature {temp}°C...")
            case_results["temperature"][f"{temp}C"] = analyze_temperature_condition(
                sequence, temp, ref_structure, ref_mea_pairs, ref_stems
            )

        # Add baseline row to parameter sweep (Turner2004)
        case_results["parameters"]["Turner2004_baseline"] = {
            "mfe_energy": ref_mfe_energy,
            "ensemble_defect": ref_ensemble_defect,
            "bp_distance_mfe": 0,
            "bp_distance_mea": 0,
            "stem_retention": baseline_retention,
        }

        # Parameter set sweep
        for param_set in ["Andronescu2007", "Langdon2018"]:
            print(f"  Parameters {param_set}...")
            case_results["parameters"][param_set] = analyze_params_condition(
                sequence, param_set, ref_structure, ref_mea_pairs, ref_stems
            )

        # Window context sweep
        print("  Window context...")
        for flank_size in [0, 25, 50, 100]:
            print(f"    Flank size {flank_size} nt...")

            if flank_size == 0:
                # No flanks: compute actual MFE-in-MEA retention (same as baseline)
                case_results["window_context"][f"flank_{flank_size}"] = {
                    "flank_size": 0,
                    "left_flank_len": 0,
                    "right_flank_len": 0,
                    "extended_len": len(sequence),
                    "mfe_energy": ref_mfe_energy,
                    "bp_distance_mea": 0,
                    "stem_retention": baseline_retention,
                    "left_accession_coords": None,
                    "right_accession_coords": None,
                    "left_clipped": False,
                    "right_clipped": False,
                }
            else:
                left_flank, right_flank = get_flanking_sequences(case_name, flank_size, cache_dir)
                if left_flank or right_flank:
                    print(f"      Fetched left={len(left_flank)} nt, right={len(right_flank)} nt")
                    case_results["window_context"][f"flank_{flank_size}"] = analyze_window_context(
                        sequence,
                        left_flank,
                        right_flank,
                        flank_size,
                        case_name,
                        ref_structure,
                        ref_mea_pairs,
                        ref_stems,
                    )
                else:
                    print("      Could not fetch flanks")

        results[case_name] = case_results

    # Save detailed JSON
    json_path = output_dir / "layer5_complete.json"
    with open(json_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n✓ Saved detailed results to {json_path}")

    # Create CSV tables
    create_temperature_tables(results, output_dir)
    create_parameters_tables(results, output_dir)
    create_window_context_tables(results, output_dir)

    return results


def create_temperature_tables(results: Dict, output_dir: Path):
    """Create CSV tables for temperature sweep."""

    # Energy table
    rows = []
    for case_name, case_data in results.items():
        row = {
            "case": case_name,
            "37C_baseline": case_data["baseline"]["mfe_energy"],
        }
        for temp_key, temp_data in case_data["temperature"].items():
            row[temp_key] = temp_data["mfe_energy"]
        rows.append(row)

    df = pd.DataFrame(rows)
    csv_path = output_dir / "temperature_mfe_energies.csv"
    df.to_csv(csv_path, index=False, float_format="%.2f")
    print(f"✓ {csv_path.name}")

    # BP distance table
    rows = []
    for case_name, case_data in results.items():
        for temp_key, temp_data in case_data["temperature"].items():
            rows.append(
                {
                    "case": case_name,
                    "temperature": temp_key,
                    "bp_distance_mfe": temp_data["bp_distance_mfe"],
                    "bp_distance_mea": temp_data["bp_distance_mea"],
                    "ensemble_defect": temp_data["ensemble_defect"],
                }
            )

    df = pd.DataFrame(rows)
    csv_path = output_dir / "temperature_bp_distances.csv"
    df.to_csv(csv_path, index=False, float_format="%.2f")
    print(f"✓ {csv_path.name}")

    # Stem retention table
    rows = []
    for case_name, case_data in results.items():
        for temp_key, temp_data in case_data["temperature"].items():
            for tier, tier_data in temp_data["stem_retention"].items():
                if tier_data["count"] > 0:
                    rows.append(
                        {
                            "case": case_name,
                            "temperature": temp_key,
                            "tier": tier.upper(),
                            "n_stems": tier_data["count"],
                            "retention": tier_data["retention"],
                        }
                    )

    df = pd.DataFrame(rows)
    csv_path = output_dir / "temperature_stem_retention.csv"
    df.to_csv(csv_path, index=False, float_format="%.4f")
    print(f"✓ {csv_path.name}")


def create_parameters_tables(results: Dict, output_dir: Path):
    """Create CSV tables for parameter set sweep."""

    # Energy table
    rows = []
    for case_name, case_data in results.items():
        row = {
            "case": case_name,
            "Turner2004_baseline": case_data["baseline"]["mfe_energy"],
        }
        for param_name, param_data in case_data["parameters"].items():
            if param_name != "Turner2004_baseline":
                row[param_name] = param_data["mfe_energy"]
        rows.append(row)

    df = pd.DataFrame(rows)
    csv_path = output_dir / "parameters_mfe_energies.csv"
    df.to_csv(csv_path, index=False, float_format="%.2f")
    print(f"✓ {csv_path.name}")

    # BP distance table
    rows = []
    for case_name, case_data in results.items():
        for param_name, param_data in case_data["parameters"].items():
            rows.append(
                {
                    "case": case_name,
                    "parameters": param_name,
                    "bp_distance_mfe": param_data["bp_distance_mfe"],
                    "bp_distance_mea": param_data["bp_distance_mea"],
                    "ensemble_defect": param_data["ensemble_defect"],
                }
            )

    df = pd.DataFrame(rows)
    csv_path = output_dir / "parameters_bp_distances.csv"
    df.to_csv(csv_path, index=False, float_format="%.2f")
    print(f"✓ {csv_path.name}")

    # Stem retention table
    rows = []
    for case_name, case_data in results.items():
        for param_name, param_data in case_data["parameters"].items():
            for tier, tier_data in param_data["stem_retention"].items():
                if tier_data["count"] > 0:
                    rows.append(
                        {
                            "case": case_name,
                            "parameters": param_name,
                            "tier": tier.upper(),
                            "n_stems": tier_data["count"],
                            "retention": tier_data["retention"],
                        }
                    )

    df = pd.DataFrame(rows)
    csv_path = output_dir / "parameters_stem_retention.csv"
    df.to_csv(csv_path, index=False, float_format="%.4f")
    print(f"✓ {csv_path.name}")


def create_window_context_tables(results: Dict, output_dir: Path):
    """Create CSV tables for window context sweep."""

    # Basic metrics table
    rows = []
    for case_name, case_data in results.items():
        if "window_context" not in case_data:
            continue

        for flank_key, flank_data in case_data["window_context"].items():
            rows.append(
                {
                    "case": case_name,
                    "flank_size": flank_data["flank_size"],
                    "left_flank_len": flank_data["left_flank_len"],
                    "right_flank_len": flank_data["right_flank_len"],
                    "extended_len": flank_data["extended_len"],
                    "mfe_energy": flank_data["mfe_energy"],
                    "bp_distance_mea": flank_data["bp_distance_mea"],
                    "left_coords": flank_data.get("left_accession_coords", ""),
                    "right_coords": flank_data.get("right_accession_coords", ""),
                    "left_clipped": flank_data.get("left_clipped", False),
                    "right_clipped": flank_data.get("right_clipped", False),
                }
            )

    df = pd.DataFrame(rows)
    csv_path = output_dir / "window_context_metrics.csv"
    df.to_csv(csv_path, index=False, float_format="%.2f")
    print(f"✓ {csv_path.name}")

    # Stem retention table
    rows = []
    for case_name, case_data in results.items():
        if "window_context" not in case_data:
            continue

        for flank_key, flank_data in case_data["window_context"].items():
            for tier, tier_data in flank_data["stem_retention"].items():
                if tier_data["count"] > 0 and tier_data["retention"] is not None:
                    rows.append(
                        {
                            "case": case_name,
                            "flank_size": flank_data["flank_size"],
                            "tier": tier.upper(),
                            "n_stems": tier_data["count"],
                            "retention": tier_data["retention"],
                        }
                    )

    df = pd.DataFrame(rows)
    csv_path = output_dir / "window_context_stem_retention.csv"
    df.to_csv(csv_path, index=False, float_format="%.4f")
    print(f"✓ {csv_path.name}")
