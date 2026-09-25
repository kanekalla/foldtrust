# Robustness Analysis

## Overview

- **N cases:** 5
- **Temperatures tested:** 24°C, 37°C, 42°C

## Baseline Tier Distribution (37°C)

| Case | Length | FIRM | SOFT | FLOPPY | Total |
|------|--------|------|------|--------|-------|
| sars2-fse | 181 | 1 | 0 | 7 | 8 |
| mapt-e10 | 268 | 7 | 10 | 2 | 19 |
| hcv-ires-dii | 268 | 12 | 4 | 2 | 18 |
| cftr-5utr | 268 | 3 | 9 | 5 | 17 |
| smn2-iss-n1 | 201 | 3 | 5 | 5 | 13 |

## Temperature Stability

Tier stability: fraction of base pairs maintaining the same FIRM/SOFT/FLOPPY classification.
Structure Jaccard: overlap of base-paired positions between baseline and test temperature.

| Case | 24°C Tier Stab. | 37°C (baseline) | 42°C Tier Stab. | 24°C Struct. Jacc. | 42°C Struct. Jacc. |
|------|-----------------|-----------------|-----------------|-------------------|-------------------|
| sars2-fse | 1.000 | 1.000 | 1.000 | 0.068 | 1.000 |
| mapt-e10 | 0.643 | 1.000 | 0.694 | 0.875 | 0.960 |
| hcv-ires-dii | 0.880 | 1.000 | 0.878 | 0.500 | 1.000 |
| cftr-5utr | 0.594 | 1.000 | 0.877 | 0.719 | 1.000 |
| smn2-iss-n1 | 0.895 | 1.000 | 0.804 | 1.000 | 0.895 |

## Summary Statistics

- **Mean FIRM stems (37°C):** 5.2
- **Mean SOFT stems (37°C):** 5.6
- **Mean FLOPPY stems (37°C):** 4.2

### Temperature Stability (mean across cases)

| Temperature | Tier Stability | Structure Jaccard |
|-------------|----------------|-------------------|
| 24C | 0.802 | 0.632 |
| 37C | 1.000 | 1.000 |
| 42C | 0.851 | 0.971 |


## Methods

- **Temperature sweep:** RNAfold with `-T` flag (24°C, 37°C, 42°C)
- **Tier stability:** Fraction of common base pairs maintaining the same tier classification
- **Structure Jaccard:** Jaccard index of base-paired positions

## Interpretation

High tier stability (close to 1.0) indicates that FIRM/SOFT/FLOPPY classifications are
robust to temperature changes. High structure Jaccard indicates that the same base pairs
form across temperatures. Disease windows with low stability may have competing folds
sensitive to physiological temperature variation.
