# Layer 3: Calibration

## Question

Are ViennaRNA base-pair probabilities well-calibrated? That is, do pairs with predicted probability *p* occur with frequency *p* in known structures?

## Data

Same 600-structure sample from Layer 2 (ArchiveII 200, bpRNA TS0 200, Rfam seed 200; all ≤500 nt, seed 42).

For each structure:
- **Reference**: Known pairs (after pseudoknot removal)
- **Candidates**: All pairs with base-pair probability > 0.001
- **MFE structure**: Used for tier analysis

Total: 181,390 candidate pairs, 20,437 correct (11.27%).

## Method

1. Fold each sequence with ViennaRNA (Turner2004, 37°C)
2. Compute partition function and base-pair probability matrix
3. For each candidate pair (i,j) with p > 0.001:
   - Record predicted probability *p*
   - Record ground truth: 1 if (i,j) is in reference, 0 otherwise
4. Compute calibration metrics:
   - **ECE** (Expected Calibration Error): mean |predicted - actual| over 10 equally-spaced bins
   - **AUROC**: Area under ROC curve (discriminative power)
   - **AUPRC**: Area under precision-recall curve
5. Tier analysis:
   - Bin MFE pairs by their probability: FIRM (p≥0.85), SOFT (0.5≤p<0.85), FLOPPY (p<0.5)
   - Report pooled PPV per tier (sum correct / sum pairs across all structures)

**Command**:
```bash
python scripts/run_layers_1_2_3.py
```

## Results

### Overall Calibration

| Metric | Value |
|--------|------:|
| ECE | 0.0692 |
| ECE (p≥0.5) | 0.3302 |
| AUROC | 0.8870 |
| AUPRC | 0.6028 |
| Candidates (p>0.001) | 181,390 |
| Positive pairs | 20,437 (11.27%) |

### Reliability Per Bin

| Bin | Mean Prob | Fraction Correct | Count |
|-----|-----------|------------------|-------|
| [0.0, 0.1) | 0.014 | 0.027 | 134,654 |
| [0.1, 0.2) | 0.143 | 0.109 | 8,937 |
| [0.2, 0.3) | 0.248 | 0.147 | 4,940 |
| [0.3, 0.4) | 0.349 | 0.188 | 3,669 |
| [0.4, 0.5) | 0.452 | 0.216 | 3,409 |
| [0.5, 0.6) | 0.548 | 0.298 | 2,732 |
| [0.6, 0.7) | 0.651 | 0.285 | 2,579 |
| [0.7, 0.8) | 0.750 | 0.344 | 2,602 |
| [0.8, 0.9) | 0.854 | 0.410 | 3,671 |
| [0.9, 1.0] | 0.976 | 0.681 | 14,197 |

### MFE-Pair Tier PPV

| Tier | Pooled PPV | N Structures | N Pairs |
|------|------------|--------------|---------|
| FIRM | 0.664 | 600 | 16,202 |
| SOFT | 0.335 | 600 | 8,191 |
| FLOPPY | 0.177 | 600 | 5,450 |

**Tier definitions**: FIRM = p≥0.85, SOFT = 0.5≤p<0.85, FLOPPY = p<0.5. MFE pairs are binned by their own base-pair probability. Pooled PPV = (sum correct pairs) / (sum all pairs) across all structures. N Structures = number of structures with ≥1 pair in that tier (same structure can contribute to multiple tiers).

## Interpretation

**Discrimination (ranking)**: ViennaRNA probabilities rank pairs well (AUROC 0.887), distinguishing correct from incorrect pairs significantly better than random (0.5).

**Calibration (accuracy)**: Probabilities are **systematically overconfident** for p ≥ 0.5:
- Pairs predicted at mean p ~ 0.95 (bin [0.9,1.0]) are correct only 68% of the time (gap: 27 points)
- Pairs predicted at p ~ 0.85 (bin [0.8,0.9)) are correct 41% of the time (gap: 44 points)
- **ECE restricted to p≥0.5 = 0.330** (overall ECE 0.069 is misleading because 74% of candidates sit in [0,0.1) where calibration is near-perfect)

**Practical impact (MFE tiers)**:
- **FIRM pairs** (mean p ~ 0.95) are correct ~66% of the time. This is **useful but far from certain**: roughly 1 in 3 high-confidence pairs is wrong.
- **SOFT pairs** (p ~ 0.5–0.85) are correct ~34% of the time
- **FLOPPY pairs** (p < 0.5) are correct ~18% of the time (slightly better than the 11% base rate, but not by much)

**Conclusion**: Base-pair probabilities successfully rank pairs, but **should not be interpreted as literal confidence**. A "95% confident" pair is actually ~66% likely to be correct. Users relying on these probabilities for decision-making (e.g., primer design, therapeutic targeting) should apply empirical recalibration.

## Limitations

1. **Reference bias**: Same as Layer 2 (comparative structures may already reflect thermodynamic assumptions)
2. **Turner2004 only**: Other parameter sets may have different calibration
3. **Sample size**: 600 structures is adequate for overall trends but not for rare-event analysis
4. **Bin discretization**: 10 bins may smooth over finer calibration patterns
5. **No recalibration attempted**: We report raw ViennaRNA probabilities; post-hoc recalibration (isotonic regression, Platt scaling) could improve calibration
6. **MFE-only tiers**: MEA and centroid tier analysis not included (future work)
