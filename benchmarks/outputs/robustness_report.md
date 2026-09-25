# Robustness Analysis

## Baseline Tier Distribution

- **N cases:** 5

| Case | Length | FIRM | SOFT | FLOPPY | Total |
|------|--------|------|------|--------|-------|
| cftr-5utr | 268 | 3 | 9 | 5 | 17 |
| hcv-ires-dii | 268 | 12 | 4 | 2 | 18 |
| mapt-e10 | 268 | 7 | 10 | 2 | 19 |
| sars2-fse | 181 | 1 | 0 | 7 | 8 |
| smn2-iss-n1 | 201 | 3 | 5 | 5 | 13 |

## Summary

- **Mean FIRM stems:** 5.2
- **Mean SOFT stems:** 5.6
- **Mean FLOPPY stems:** 4.2

## Limitations

Temperature and parameter set sweep requires ViennaRNA Python library (future enhancement)

For full robustness testing (temperature sweep, parameter sets, window jitter), 
the ViennaRNA Python library would be needed to programmatically control fold 
parameters. The CLI-based approach used here is limited to default parameters.

**Future work:**
- Temperature sweep (24°C, 37°C, 42°C)
- Parameter sets (Turner 2004, Andronescu 2007, Langdon 2018)
- Window boundary jitter (+/- 10-25 nt)
