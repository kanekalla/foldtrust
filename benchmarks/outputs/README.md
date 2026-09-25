# FoldTrust Benchmark Outputs

This directory contains results from Layers 1-3 of the FoldTrust benchmark.

## Structure

```
outputs/
├── layer1/                          # Layer 1: Scoring correctness tests
│   └── layer1_tests.json            # Parser and tier regression test results
├── layer2/                          # Layer 2: Structure accuracy
│   ├── layer2_raw_results.csv       # Per-structure predictions (all methods)
│   ├── layer2_summary_dataset.csv   # F1 by dataset and method
│   ├── layer2_summary_family.csv    # F1 by Rfam family
│   ├── layer2_summary_length.csv    # F1 by length bin
│   ├── layer2_summary_overall.csv   # Overall metrics per method
│   └── layer2_summary.json          # JSON summary with headline numbers
├── layer3/                          # Layer 3: Calibration and tier analysis
│   ├── layer3_calibration_curve.csv # Reliability diagram data
│   ├── layer3_tier_analysis.csv     # Per-structure tier PPVs
│   ├── layer3_tier_summary.csv      # Aggregated tier PPVs with CIs
│   └── layer3_summary.json          # Overall calibration metrics
├── figures/                         # Generated plots
│   ├── layer2_method_comparison.png
│   ├── layer2_f1_vs_length.png
│   ├── layer2_f1_by_family.png
│   ├── layer3_reliability_diagram.png
│   └── layer3_tier_ppv.png
└── benchmark_run.log                # Full benchmark output log
```

## Running the Benchmark

### Default (Fast)

Max length 500 nt, 200 structures per dataset, typical runtime ~45 minutes on a 16GB machine:

```bash
python3 scripts/run_layers_1_2_3.py
```

### Full (Slow)

All sequences, no length cap:

```bash
python3 scripts/run_layers_1_2_3.py --full
```

## Exact Commands Run

The benchmark runner orchestrates all three layers:

```bash
python3 scripts/run_layers_1_2_3.py
```

This calls:
- `foldtrust.benchmark.layer1_scoring.run_layer1_tests()` → `benchmarks/outputs/layer1/`
- `foldtrust.benchmark.layer2_accuracy.run_layer2_benchmark()` → `benchmarks/outputs/layer2/`
- `foldtrust.benchmark.layer3_calibration.run_layer3_calibration()` → `benchmarks/outputs/layer3/`

## Data Sources

All data is cached in `data/_cache/` (gitignored).

- **ArchiveII**: Sloma & Mathews (2016), RNA 22:1808-1818
- **bpRNA TS0**: Singh et al. (2019), Nat Commun 10:5407 via MXfold2 Zenodo
- **Rfam seed**: Ontiveros-Palacios et al. (2025), NAR 53(D1):D258-D267
- **SPOT-RNA PDB**: Singh et al. (2019), Nat Commun 10:5407 (not yet integrated)

See `data/_cache/MANIFEST.md` for full provenance.

## Documentation

Detailed documentation:
- `docs/benchmark/layer1_scoring.md` - Scoring correctness tests
- `docs/benchmark/layer2_accuracy.md` - Structure accuracy evaluation
- `docs/benchmark/layer3_calibration.md` - Calibration and tier analysis

## Headline Metrics

Generated from benchmark run on 2026-09-25 (ViennaRNA 2.7.2, Turner2004, 37°C, 600 structures).

### Layer 2: Structure Accuracy

- Total structures: 600 (ArchiveII 200, bpRNA TS0 200, Rfam seed 200)
- MFE F1: 0.548 [95% CI: 0.526, 0.570]
- MEA F1: 0.563 [95% CI: 0.542, 0.585]
- Centroid F1: 0.570 [95% CI: 0.549, 0.590]

### Layer 3: Calibration

- Candidates (p>0.001): 181,390 pairs (20,437 correct, 11.27%)
- ECE: 0.069 (overall), 0.330 (restricted to p≥0.5)
- AUROC: 0.887
- AUPRC: 0.603
- FIRM PPV (p≥0.85): 0.664 (16,202 pairs)
- SOFT PPV (0.5≤p<0.85): 0.335 (8,191 pairs)
- FLOPPY PPV (p<0.5): 0.177 (5,450 pairs)

## Integrity

All numbers in this directory come from runs performed during benchmark execution. No results are hand-edited or cherry-picked. References are never modified toward predictions.

See `data/_cache/SHA256SUMS` for data integrity verification.
