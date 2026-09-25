# Reference Structure Accuracy Benchmark

## Dataset
- **Source:** RNA STRAND / ArchiveII subset
- **N sequences:** 15
- **Allow slip:** True

## Overall Results (MFE Structure)

| Metric | Mean | Std |
|--------|------|-----|
| Sensitivity | 0.213 | 0.410 |
| PPV | 0.213 | 0.410 |
| F1 | 0.213 | 0.410 |
| MCC | 0.206 | 0.414 |

## Stratified by Length

| Length Bin | N | Sensitivity | PPV | F1 | MCC |
|------------|---|-------------|-----|----|----|
| 0-50 | 6 | 0.333 | 0.333 | 0.333 | 0.323 |
| 50-100 | 6 | 0.200 | 0.199 | 0.199 | 0.193 |
| 100-200 | 3 | 0.000 | 0.000 | 0.000 | -0.005 |

## Methods

- **MFE prediction:** ViennaRNA `RNAfold` (Turner 2004 parameters, 37°C)
- **Metrics:**
  - **Sensitivity (Recall):** TP / (TP + FN)
  - **PPV (Precision):** TP / (TP + FP)
  - **F1:** 2 * (Sensitivity * PPV) / (Sensitivity + PPV)
  - **MCC:** Matthews Correlation Coefficient
  - **Slip tolerance:** True (1-nt slippage allowed)

## Interpretation

These metrics evaluate how well ViennaRNA's MFE structure predicts the reference 
(comparative or experimentally determined) secondary structure. F1 and MCC provide 
balanced measures of prediction quality.
