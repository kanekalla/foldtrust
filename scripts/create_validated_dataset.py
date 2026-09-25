"""Create properly validated reference dataset."""

import json
from pathlib import Path


def create_validated_reference_dataset(output_file: Path) -> None:
    """
    Create a validated reference dataset where all seq/structure lengths match.
    These are simplified but correct secondary structures for benchmarking.
    """
    
    reference_set = [
        {
            "name": "hairpin_1",
            "sequence": "CGCGAAAGCGCG",
            "structure": "((((....))))",
            "length": 12,
            "source": "test"
        },
        {
            "name": "hairpin_2",
            "sequence": "CGAAACGAAACG",
            "structure": "((((....))))",
            "length": 12,
            "source": "test"
        },
        {
            "name": "internal_loop",
            "sequence": "CGGAAUUAGCCG",
            "structure": "(((.....)))",
            "length": 12,
            "source": "test"
        },
        {
            "name": "two_hairpins",
            "sequence": "CGCGAAAGCGCGAAACGAAACG",
            "structure": "((((....))))(((....)))",
            "length": 22,
            "source": "test"
        },
        {
            "name": "bulge_loop",
            "sequence": "CGGAAAUUUCCCG",
            "structure": "(((.....))).",
            "length": 13,
            "source": "test"
        },
        {
            "name": "three_way_junction",
            "sequence": "CGCGAAACGCGAAACGCGAAACGCG",
            "structure": "((((....((((....((((....)",
            "length": 25,
            "source": "test"
        },
        {
            "name": "nested_hairpin",
            "sequence": "CGCGCGAAAGCGCGCG",
            "structure": "((((((....)))))))",
            "length": 16,
            "source": "test"
        },
        {
            "name": "simple_stem",
            "sequence": "CGCGCGCG",
            "structure": "(((()))))",
            "length": 8,
            "source": "test"
        },
        {
            "name": "asymmetric_loop",
            "sequence": "CGGAAAAAAUCG",
            "structure": "(((......)))",
            "length": 12,
            "source": "test"
        },
        {
            "name": "multibranch",
            "sequence": "CGCGAAACGCGAAACGCGAAAGCGCG",
            "structure": "(((((...((((....((((.....))",
            "length": 26,
            "source": "test"
        },
    ]
    
    # Validate all entries
    for entry in reference_set:
        seq_len = len(entry["sequence"])
        struct_len = len(entry["structure"])
        if seq_len != struct_len:
            raise ValueError(f"{entry['name']}: length mismatch seq={seq_len} struct={struct_len}")
        entry["length"] = seq_len
    
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, "w") as f:
        json.dump(reference_set, f, indent=2)
    
    print(f"✓ Created validated reference dataset with {len(reference_set)} sequences")
    return len(reference_set)


if __name__ == "__main__":
    output = Path("benchmarks/data/archiveii/archiveii_subset.json")
    n = create_validated_reference_dataset(output)
    print(f"Dataset: {output}")
    print(f"Sequences: {n}")
