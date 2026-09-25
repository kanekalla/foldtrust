# Layer 3: Calibration

## Question

Are ViennaRNA base-pair probabilities well-calibrated? That is, do pairs with predicted probability *p* occur with frequency *p* in known structures?

## Data

Same 600-structure sample from Layer 2 (ArchiveII 200, bpRNA TS0 200, Rfam seed 200; all ≤500 nt, seed 42).

For each structure:
- **Reference**: Known pairs (after greedy pseudoknot removal, same as Layer 2)
- **Candidates**: All pairs (i<j) with base-pair probability p > 0.001
- **MFE structure**: Used for tier analysis

**Definitions**:
- **Candidates**: All pairs i<j with bpp p > 1e-3
- **Positives**: Pairs in the reference after the same greedy pseudoknot removal as Layer 2
- **10 equal-width bins** on [0,1] (last bin closed)
- **ECE** = sum_b (n_b/N) |mean_p_b - frac_true_b|
- **ECE restricted to p≥0.5** = same formula over only candidates with p≥0.5
- **AUROC** over all candidates
- **Tier PPV over MFE pairs**: MFE pairs binned by their bpp into FIRM p≥0.85 / SOFT 0.5-0.85 / FLOPPY <0.5; pooled PPV = true/total

Total: 185,653 candidate pairs, 21,191 correct (11.41%).

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
| ECE | 0.0664 |
| ECE (p≥0.5) | 0.3174 |
| AUROC | 0.8889 |
| AUPRC | 0.6153 |
| Candidates (p>0.001) | 185,653 |
| Positive pairs | 21,191 (11.41%) |

### Reliability Per Bin

| Bin | Mean Prob | Fraction Correct | Count |
|-----|-----------|------------------|-------|
| [0.0, 0.1) | 0.014 | 0.027 | 138,095 |
| [0.1, 0.2) | 0.143 | 0.108 | 9,138 |
| [0.2, 0.3) | 0.248 | 0.156 | 5,027 |
| [0.3, 0.4) | 0.350 | 0.187 | 3,698 |
| [0.4, 0.5) | 0.452 | 0.224 | 3,359 |
| [0.5, 0.6) | 0.548 | 0.279 | 2,695 |
| [0.6, 0.7) | 0.650 | 0.278 | 2,611 |
| [0.7, 0.8) | 0.750 | 0.358 | 2,839 |
| [0.8, 0.9) | 0.853 | 0.441 | 3,833 |
| [0.9, 1.0] | 0.976 | 0.700 | 14,358 |

### MFE-Pair Tier PPV

| Tier | Pooled PPV | N Structures | N Pairs |
|------|------------|--------------|---------|
| FIRM | 0.674 | 582 | 16,417 |
| SOFT | 0.299 | 585 | 8,518 |
| FLOPPY | 0.146 | 466 | 5,504 |

**Tier definitions**: FIRM = p≥0.85, SOFT = 0.5≤p<0.85, FLOPPY = p<0.5. MFE pairs are binned by their own base-pair probability. Pooled PPV = (sum correct pairs) / (sum all pairs) across all structures. N Structures = number of structures with ≥1 pair in that tier (same structure can contribute to multiple tiers).

## Interpretation

**Discrimination (ranking)**: ViennaRNA probabilities rank pairs well (AUROC 0.887), distinguishing correct from incorrect pairs significantly better than random (0.5).

**Calibration (accuracy)**: Probabilities are **systematically overconfident** for p ≥ 0.5:
- The [0.9,1.0] bin (mean p ~ 0.98) has frac_true ~0.70, i.e. high-probability pairs are overconfident against these references
- Pairs predicted at p ~ 0.85 (bin [0.8,0.9)) are correct 44% of the time (gap: 41 points)
- **ECE restricted to p≥0.5 = 0.317** (overall ECE 0.066 is dominated by the large [0,0.1) bin where calibration is near-perfect)

**Practical impact (MFE tiers)**:
- **FIRM pairs** (p≥0.85) are correct ~67% of the time. This is **useful but far from certain**: roughly 1 in 3 high-confidence pairs is wrong.
- **SOFT pairs** (p ~ 0.5–0.85) are correct ~30% of the time
- **FLOPPY pairs** (p < 0.5, MFE only) are correct ~15% of the time (vs 4.4% for all-candidate FLOPPY)

**Conclusion**: Base-pair probabilities successfully rank pairs, but **should not be interpreted as literal confidence**. A "95% confident" pair is actually ~70% likely to be correct against these references. Users relying on these probabilities for decision-making (e.g., primer design, therapeutic targeting) should apply empirical recalibration.

## Limitations

1. **Reference bias**: Same as Layer 2 (comparative structures may already reflect thermodynamic assumptions)
2. **Turner2004 only**: Other parameter sets may have different calibration
3. **Sample size**: 600 structures is adequate for overall trends but not for rare-event analysis
4. **Bin discretization**: 10 bins may smooth over finer calibration patterns
5. **No recalibration attempted**: We report raw ViennaRNA probabilities; post-hoc recalibration (isotonic regression, Platt scaling) could improve calibration
6. **MFE-only tiers**: MEA and centroid tier analysis not included (future work)

