#!/usr/bin/env python3
"""Extract case metrics from FoldTrust analysis for Layer 0 documentation.

This script:
1. Runs FoldTrust analysis on each case
2. Extracts key metrics (MFE, ensemble free energy, stem counts, GC%)
3. Saves to CSV and JSON for documentation
"""

import json
import sys
from pathlib import Path

import numpy as np
import yaml

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

import RNA

from foldtrust import vienna


def compute_gc_content(sequence: str) -> float:
    """Compute GC content as percentage."""
    gc_count = sequence.upper().count("G") + sequence.upper().count("C")
    return 100.0 * gc_count / len(sequence)


def count_stems_by_tier(prob_matrix: np.ndarray, structure: str, min_stem_length: int = 3) -> dict:
    """Count stems by reliability tier (FIRM >= 0.85, SOFT 0.5-0.85, FLOPPY < 0.5)."""
    # Parse pairs from structure
    pairs = []
    stack = []

    for i, char in enumerate(structure):
        if char in "([{":
            stack.append(i)
        elif char in ")]}":
            if stack:
                j = stack.pop()
                pairs.append((j, i))

    pairs.sort()

    # Group consecutive base pairs into stems
    stems = []
    current_stem = []

    for i, j in pairs:
        if current_stem and (i != current_stem[-1][0] + 1 or j != current_stem[-1][1] - 1):
            # Start new stem
            if len(current_stem) >= min_stem_length:
                stems.append(current_stem)
            current_stem = [(i, j)]
        else:
            current_stem.append((i, j))

    if len(current_stem) >= min_stem_length:
        stems.append(current_stem)

    # Classify stems by minimum pair probability
    firm_count = 0
    soft_count = 0
    floppy_count = 0

    for stem in stems:
        min_prob = min(prob_matrix[i, j] for i, j in stem)

        if min_prob >= 0.85:
            firm_count += 1
        elif min_prob >= 0.5:
            soft_count += 1
        else:
            floppy_count += 1

    return {"firm": firm_count, "soft": soft_count, "floppy": floppy_count, "total": len(stems)}


def get_verdict(stem_counts: dict) -> str:
    """Determine verdict from stem counts (simplified logic)."""
    if stem_counts["floppy"] > 0:
        return "REDESIGN"
    elif stem_counts["soft"] > stem_counts["firm"]:
        return "NEED PROBING"
    else:
        return "TRUST"


def analyze_case(case_dir: Path) -> dict:
    """Analyze one case and return metrics."""
    case_name = case_dir.name

    sequence_file = case_dir / "sequence.fa"
    meta_file = case_dir / "meta.yaml"

    # Read sequence
    with open(sequence_file) as f:
        lines = f.readlines()
    sequence = "".join(line.strip() for line in lines[1:] if not line.startswith(">"))

    # Read metadata
    meta = yaml.safe_load(meta_file.read_text())
    prov = meta["provenance"]

    # Compute structure (Turner 2004, default 37°C)
    RNA.params_load_RNA_Turner2004()
    md = RNA.md()
    md.temperature = 37.0

    structure, mfe = vienna.fold_mfe(sequence)

    # Compute partition function
    fc = RNA.fold_compound(sequence, md)
    ensemble_energy = fc.pf()[1]  # Returns (structure, free_energy)

    # Compute pair probabilities
    pair_probs = vienna.compute_pair_probabilities(sequence)

    # Compute metrics
    gc_content = compute_gc_content(sequence)
    stem_counts = count_stems_by_tier(pair_probs, structure)
    verdict = get_verdict(stem_counts)

    # Extract coordinates
    coord_str = f"{prov['accession']}:{prov['start']}-{prov['end']}({prov['strand']})"

    return {
        "case": case_name,
        "coordinates": coord_str,
        "accession": prov["accession"],
        "start": prov["start"],
        "end": prov["end"],
        "strand": prov["strand"],
        "length": len(sequence),
        "gc_percent": round(gc_content, 1),
        "mfe_kcal_mol": round(mfe, 2),
        "mfe_structure": structure,
        "ensemble_free_energy": round(ensemble_energy, 2),
        "firm_stems": stem_counts["firm"],
        "soft_stems": stem_counts["soft"],
        "floppy_stems": stem_counts["floppy"],
        "total_stems": stem_counts["total"],
        "verdict": verdict,
    }


def main():
    cases_dir = Path("data/cases")
    output_dir = Path("benchmarks/outputs/layer0_cases")
    output_dir.mkdir(parents=True, exist_ok=True)

    case_dirs = sorted([d for d in cases_dir.iterdir() if d.is_dir()])

    results = []

    print("Analyzing cases...")
    for case_dir in case_dirs:
        print(f"  {case_dir.name}...", end=" ", flush=True)
        metrics = analyze_case(case_dir)
        results.append(metrics)
        print(f"MFE={metrics['mfe_kcal_mol']:.2f}, verdict={metrics['verdict']}")

    # Save CSV
    csv_file = output_dir / "case_metrics.csv"
    with open(csv_file, "w") as f:
        # Header
        headers = [
            "case",
            "coordinates",
            "length",
            "gc_percent",
            "mfe_kcal_mol",
            "mfe_structure",
            "ensemble_free_energy",
            "firm_stems",
            "soft_stems",
            "floppy_stems",
            "total_stems",
            "verdict",
        ]
        f.write(",".join(headers) + "\n")

        # Data
        for r in results:
            row = [str(r[h]) for h in headers]
            f.write(",".join(row) + "\n")

    print(f"\nSaved CSV to {csv_file}")

    # Save JSON
    json_file = output_dir / "case_metrics.json"
    with open(json_file, "w") as f:
        json.dump(results, f, indent=2)

    print(f"Saved JSON to {json_file}")


if __name__ == "__main__":
    main()
