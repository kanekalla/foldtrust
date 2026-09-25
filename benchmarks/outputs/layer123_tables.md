# Layer 1-3 Benchmark Tables

Generated from benchmark outputs in `benchmarks/outputs/layer2/` and `layer3/`.

## Layer 2: Structure Accuracy

### Overall Metrics (Mean of Per-Structure)

| Method | Sensitivity | PPV | F1 | MCC |
|--------|-------------|-----|----|----|
| MFE | 0.639 | 0.499 | 0.548 | 0.556 |
| MEA | 0.641 | 0.521 | 0.563 | 0.569 |
| CENTROID | 0.622 | 0.551 | 0.570 | 0.576 |

### Per-Dataset F1 (Mean of Per-Structure)

| Dataset | N | MFE | MEA | Centroid |
|---------|---|----|-----|---------
| ArchiveII | 200 | 0.562 | 0.582 | 0.591 |
| Rfam_seed | 200 | 0.569 | 0.582 | 0.586 |
| bpRNA_TS0 | 200 | 0.514 | 0.525 | 0.532 |

### Sample Description

- **Total structures**: 600

- **Sampling**: Random seed 42, 200 structures per dataset (after ≤500 nt filter), re-seeded before each dataset

- **Pool sizes**: ArchiveII 3854, bpRNA TS0 1305, Rfam seed 468 (all ≤500 nt)

- **Duplicates**: Not removed (9 duplicate-sequence pairs exist in the 600)

- **Pseudoknot handling**: Greedy nested-only selection removes crossing pairs from reference

- **Matching**: Exact base-pair matching (i,j exact)

- **Slip tolerance**: Predicted pair (i,j) matches if any of (i±1,j), (i,j±1), or (i,j) is in reference; 
reference pair recovered if any of (i±1,j), (i,j±1), or (i,j) is predicted. 
Slip F1 ≥ exact F1 by construction.

- **Averaging**: F1 = mean of per-structure F1 (primary); bootstrap 95% CI with 1000 resamples


## Layer 3: Calibration

### Overall Calibration

| Metric | Value |
|--------|------:|
| ECE | 0.0692 |
| ECE (p≥0.5) | 0.3302 |
| AUROC | 0.8870 |
| AUPRC | 0.6028 |
| Candidates (p>0.001) | 181390 |
| Positive pairs | 20437 (11.27%) |

### Reliability Per Bin

| Bin | Mean Prob | Fraction Correct | Count |
|-----|-----------|------------------|-------|
| [0.0, 0.1) | 0.014 | 0.027 | 134654 |
| [0.1, 0.2) | 0.143 | 0.109 | 8937 |
| [0.2, 0.3) | 0.248 | 0.147 | 4940 |
| [0.3, 0.4) | 0.349 | 0.188 | 3669 |
| [0.4, 0.5) | 0.452 | 0.216 | 3409 |
| [0.5, 0.6) | 0.548 | 0.298 | 2732 |
| [0.6, 0.7) | 0.651 | 0.285 | 2579 |
| [0.7, 0.8) | 0.750 | 0.344 | 2602 |
| [0.8, 0.9) | 0.854 | 0.410 | 3671 |
| [0.9, 1.0] | 0.976 | 0.681 | 14197 |

### MFE-Pair Tier PPV

| Tier | Pooled PPV | N Structures | N Pairs |
|------|------------|--------------|---------|
| FIRM | 0.664 | 600 | 16202 |
| SOFT | 0.335 | 600 | 8191 |
| FLOPPY | 0.177 | 600 | 5450 |

**Tier definitions**: FIRM = p≥0.85, SOFT = 0.5≤p<0.85, FLOPPY = p<0.5. 
MFE pairs are binned by their own base-pair probability. 
Pooled PPV = (sum correct pairs) / (sum all pairs) across all structures in tier. 
N Structures = number of structures with ≥1 pair in that tier (note: same structure can contribute to multiple tiers).