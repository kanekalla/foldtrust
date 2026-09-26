# Layer 3: Calibration

Calibration of ViennaRNA base-pair probabilities: do pairs with probability *p* occur with frequency *p*?

## Command

```bash
python3 scripts/run_layers_1_2_3.py
```

This runs `foldtrust.benchmark.layer3_calibration.run_layer3_calibration()`.

## Data

Same 600-structure sample from Layer 2. All candidate pairs with p > 0.001 are evaluated against reference structures.

Total: 181,390 candidates, 20,437 correct (11.27%).

## Outputs

- `layer3_summary.json`: Overall metrics (ECE, AUROC, AUPRC, tier PPVs with CIs)
- `layer3_calibration_curve.csv`: Reliability diagram data (10 bins)
- `layer3_tier_summary.csv`: Aggregated tier PPVs (FIRM/SOFT/FLOPPY) with bootstrap CIs
- `layer3_tier_analysis.csv`: Per-structure tier PPVs

## Headline Results

- ECE: 0.069 (overall), 0.330 (p≥0.5) — **overconfident for high probabilities**
- AUROC: 0.887 (good discrimination)
- FIRM PPV (p≥0.85): 0.674 (16,417 MFE pairs; layer3_mfe_tier_summary.csv)
- SOFT PPV (0.5≤p<0.85): 0.299 (8,518 pairs)
- FLOPPY PPV (p<0.5): 0.146 (5,504 pairs)
- Legacy files: layer3_tier_summary.csv and layer3_tier_analysis.csv (from earlier pipeline version)

**Key finding**: Probabilities rank pairs well but are systematically overconfident. A "95% confident" pair is actually ~67% correct.

See `docs/benchmark/layer3_calibration.md` for full results and interpretation.
