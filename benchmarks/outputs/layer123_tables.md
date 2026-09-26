# Layer 1-3 Benchmark Tables

Generated from benchmark outputs in `benchmarks/outputs/layer2/` and `layer3/`.

## Layer 2: Structure Accuracy

### Overall Metrics (Mean of Per-Structure)

| Method | Sensitivity | PPV | F1 | MCC |
|--------|-------------|-----|----|----|
| MFE | 0.639 | 0.499 | 0.548 | 0.557 |
| MEA | 0.641 | 0.521 | 0.563 | 0.571 |
| CENTROID | 0.622 | 0.551 | 0.570 | 0.577 |

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
| ECE | 0.0664 |
| ECE (p≥0.5) | 0.3174 |
| AUROC | 0.8889 |
| AUPRC | 0.6153 |
| Candidates (p>0.001) | 185653 |
| Positive pairs | 21191 (11.41%) |

### Reliability Per Bin

| Bin | Mean Prob | Fraction Correct | Count |
|-----|-----------|------------------|-------|
| [0.0, 0.1) | 0.014 | 0.027 | 138095 |
| [0.1, 0.2) | 0.143 | 0.108 | 9138 |
| [0.2, 0.3) | 0.248 | 0.156 | 5027 |
| [0.3, 0.4) | 0.350 | 0.187 | 3698 |
| [0.4, 0.5) | 0.452 | 0.224 | 3359 |
| [0.5, 0.6) | 0.548 | 0.279 | 2695 |
| [0.6, 0.7) | 0.650 | 0.278 | 2611 |
| [0.7, 0.8) | 0.750 | 0.358 | 2839 |
| [0.8, 0.9) | 0.853 | 0.441 | 3833 |
| [0.9, 1.0] | 0.976 | 0.700 | 14358 |

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