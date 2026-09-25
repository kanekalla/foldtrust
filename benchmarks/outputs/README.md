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

Max length 500 nt, typical runtime ~10-20 minutes:

```bash
cd /workspace
python3 scripts/run_layers_1_2_3.py
```

### Full (Slow)

All sequences, no length cap:

```bash
python3 scripts/run_layers_1_2_3.py --full
```

## Exact Commands Run

Layer 1:
```bash
python3 src/foldtrust/benchmark/layer1_scoring.py
```

Layer 2:
```bash
# Loads ArchiveII, bpRNA TS0, Rfam seed
# Evaluates MFE, MEA, centroid
# Outputs per-structure results and summaries
python3 src/foldtrust/benchmark/layer2_accuracy.py
```

Layer 3:
```bash
# Pools candidate pairs from Layer 2 datasets
# Computes calibration metrics and tier PPVs
python3 src/foldtrust/benchmark/layer3_calibration.py
```

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

*To be populated after benchmark completes.*

### Layer 2: Structure Accuracy

- Total structures: [TBD]
- MFE F1: [TBD] [95% CI]
- MEA F1: [TBD] [95% CI]

### Layer 3: Calibration

- ECE: [TBD]
- AUROC: [TBD]
- FIRM PPV: [TBD] [95% CI]
- SOFT PPV: [TBD] [95% CI]
- FLOPPY PPV: [TBD] [95% CI]

## Integrity

All numbers in this directory come from runs performed during benchmark execution. No results are hand-edited or cherry-picked. References are never modified toward predictions.

See `data/_cache/SHA256SUMS` for data integrity verification.
