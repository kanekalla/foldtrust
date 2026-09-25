# Layer 2: Structure Accuracy

Accuracy of ViennaRNA MFE, MEA, and centroid structures vs. known reference structures.

## Command

```bash
python3 scripts/run_layers_1_2_3.py
```

This runs `foldtrust.benchmark.layer2_accuracy.run_layer2_benchmark()`.

## Datasets

- ArchiveII (Mathews lab): 3,854 structures ≤500 nt → 200 sampled
- bpRNA TS0: 1,305 structures ≤500 nt → 200 sampled
- Rfam seed (15.1): 468 structures ≤500 nt → 200 sampled

Sampling: seed 42, re-seeded before each dataset.

## Outputs

- `sample_ids.csv`: List of 600 sampled structures (dataset, name, length, seq_sha1)
- `layer2_raw_results.csv`: Per-structure predictions and metrics (all methods)
- `layer2_summary_overall.csv`: Overall metrics (sensitivity, PPV, F1, MCC) with 95% CI
- `layer2_summary_dataset.csv`: F1 by dataset and method
- `layer2_summary_family.csv`: F1 by Rfam family (for Rfam subset)
- `layer2_summary_length.csv`: F1 by length bin
- `layer2_summary.json`: JSON with headline numbers

## Headline Results

- MFE F1: 0.548
- MEA F1: 0.563
- Centroid F1: 0.570 (best)

See `docs/benchmark/layer2_accuracy.md` for full results and interpretation.
