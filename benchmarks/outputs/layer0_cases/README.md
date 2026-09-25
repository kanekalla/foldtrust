# Layer 0: Case Metrics

This directory contains metrics computed from FoldTrust analysis of the five disease-relevant RNA windows.

## Files

- `case_metrics.csv` - Case-by-case metrics in CSV format
- `case_metrics.json` - Same data in JSON format

## Command

Metrics computed with:

```bash
python3 scripts/compute_case_metrics.py
```

## Software Versions

- ViennaRNA: 2.7.2 (Python package, via `pip install ViennaRNA`)
- Parameters: Turner 2004 (default)
- Temperature: 37°C (default)

## Methodology

For each case:

1. Load sequence from `data/cases/<case>/sequence.fa`
2. Compute MFE structure with `RNA.fold_compound(seq).mfe()`
3. Compute ensemble free energy with `RNA.fold_compound(seq).pf()`
4. Compute base-pair probability matrix
5. Parse stems from MFE structure (min 3 consecutive base pairs)
6. Classify stems by minimum pair probability:
   - **FIRM**: min prob ≥ 0.85 (high ensemble support)
   - **SOFT**: 0.5 ≤ min prob < 0.85 (moderate support)
   - **FLOPPY**: min prob < 0.5 (low support, fluctuates in ensemble)
7. Assign verdict:
   - **REDESIGN**: contains floppy stems
   - **NEED PROBING**: soft stems outnumber firm
   - **TRUST**: strong ensemble support

## Notes

All MFE values match the expected sanity checks specified in the task requirements:
- smn2-iss-n1: -31.30 kcal/mol ✓
- cftr-5utr: -60.60 kcal/mol ✓
- mapt-e10: -52.30 kcal/mol ✓
- hcv-ires-dii: -23.80 kcal/mol ✓
- sars2-fse: -26.00 kcal/mol ✓

These values are computed using Turner 2004 parameters at 37°C via the ViennaRNA Python API (2.7.2).
