# Reference Structure Accuracy Benchmark

## Dataset
- **Source:** RNA STRAND / ArchiveII subset
- **N sequences:** 15
- **Allow slip:** False

## Overall Results (MFE Structure)

| Metric | Mean | Std |
|--------|------|-----|
| Sensitivity | 0.067 | 0.258 |
| PPV | 0.067 | 0.258 |
| F1 | 0.067 | 0.258 |
| MCC | 0.047 | 0.264 |

## Stratified by Length

| Length Bin | N | Sensitivity | PPV | F1 | MCC |
|------------|---|-------------|-----|----|----|
| 0-50 | 12 | 0.083 | 0.083 | 0.083 | 0.062 |
| 50-100 | 3 | 0.000 | 0.000 | 0.000 | -0.010 |

## Methods

- **MFE prediction:** ViennaRNA `RNAfold` (Turner 2004 parameters, 37°C)
- **Metrics:**
  - **Sensitivity (Recall):** TP / (TP + FN)
  - **PPV (Precision):** TP / (TP + FP)
  - **F1:** 2 * (Sensitivity * PPV) / (Sensitivity + PPV)
  - **MCC:** Matthews Correlation Coefficient
  - **Slip tolerance:** False (1-nt slippage allowed)

## Interpretation

These metrics evaluate how well ViennaRNA's MFE structure predicts the reference 
(comparative or experimentally determined) secondary structure. F1 and MCC provide 
balanced measures of prediction quality.
