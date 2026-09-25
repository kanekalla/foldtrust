# Calibration Analysis

## Base-Pair Probability Calibration

- **N sequences:** 10

| Metric | Mean | Std |
|--------|------|-----|
| Expected Calibration Error (ECE) | 0.0336 | 0.0136 |
| AUROC | 0.521 | 0.107 |
| AUPRC | 0.049 | 0.082 |

**Interpretation:**
- ECE measures the average gap between predicted probabilities and observed accuracy.
- AUROC and AUPRC measure how well pair probabilities discriminate true from false pairs.

## Reliability Tier Accuracy

| Tier | Mean PPV | Total Pairs |
|------|----------|-------------|
| **FIRM** (P ≥ 0.85) | 0.000 | 45 |
| **SOFT** (0.5 ≤ P < 0.85) | 0.000 | 34 |
| **FLOPPY** (P < 0.5) | 0.000 | 16 |

**Interpretation:**
- PPV (Positive Predictive Value) shows the fraction of pairs in each tier that match the reference.
- Higher PPV for FIRM tier validates FoldTrust's reliability labeling.
- Lower PPV for FLOPPY tier confirms that low-probability pairs are less trustworthy.

## Methods

- **Calibration:** Base-pair probabilities binned into 10 intervals; observed accuracy computed per bin.
- **ECE:** Sum of |predicted - observed| weighted by bin frequency.
- **Tier classification:** Stems classified by mean pair probability (FIRM/SOFT/FLOPPY).

See `reliability_diagram.png` for visual calibration plot.
