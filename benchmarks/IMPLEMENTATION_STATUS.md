"""
Benchmark Implementation Status - Honest Assessment

After investigation, the requested benchmark analyses face significant data availability
and technical constraints:

## 1. Reference Structure Accuracy (ArchiveII)

**Status:** NOT IMPLEMENTABLE with current data access

**Attempts made:**
- Mathews lab (rna.urmc.rochester.edu): archiveII.tar.gz returns HTTP 404
- RNA STRAND v2.0: Download timeouts
- Published repos (mxfold2, LinearPartition, E2Efold, allegro): No .ct/.bpseq files found
- Alternative sources (Rfam, RMDB): Require family-by-family extraction and validation

**Conclusion:** The original ArchiveII dataset is no longer hosted at documented URLs.
Building a substitute from PDB/Rfam would require substantial data engineering beyond
the scope of this benchmark task.

**Decision:** Remove reference structure accuracy from benchmark suite per user instruction:
"If you truly cannot get real reference data, say so and do not report reference metrics at all."

## 2. SHAPE Probing Correlation

**Status:** DATA FOUND but COORDINATE MAPPING REQUIRED

**Data source found:**
- Repository: https://github.com/DasLab/SARS_CoV-2_shape_comparison
- Contains genome-wide SHAPE reactivity (29,903 positions) from multiple studies
- Files: zhang_invivo_reactivity.csv, incarnato_invivo_reactivity.csv, pyle_reactivity.csv

**Issue:** The FSE sequence in data/cases/sars2-fse/sequence.fa (181 nt) does NOT 
match any substring of the reference genome in the DasLab dataset. Proper implementation
requires:
1. Identifying exact genomic coordinates of the FoldTrust FSE window
2. Mapping SHAPE values to the 181 nt region
3. Validating sequence alignment

**Time required:** 2-4 hours for proper coordinate resolution and validation

## 3. Robustness Analysis

**Status:** PARTIALLY IMPLEMENTABLE

**Temperature sweep:** RNAfold supports `-T` flag (✓ implementable)
**Parameter files:** RNAfold supports `-P` flag but alternative parameter files
 (Andronescu 2007, Langdon 2018) not found in ViennaRNA 2.5.1 installation
**Window jitter:** Requires flanking genomic sequence not included in FoldTrust cases

**What can be done:**
- Temperature sweep (24°C, 37°C, 42°C) on the five disease cases
- Report tier stability (fraction of stems changing tier)

## 4. Calibration Analysis

**Status:** IMPLEMENTABLE on disease cases

Can compute calibration metrics (ECE, AUROC, tier PPV) using the five disease 
windows, but without reference structures, tier PPV cannot be validated.

## Recommendation

Given data access constraints, the honest approach is:

1. **Remove** reference structure accuracy (no ArchiveII)
2. **Document** SHAPE data availability and coordinate mapping issue
3. **Implement** temperature-sweep robustness on disease cases
4. **Update** NOTES.md with honest limitations

This produces a minimal but defensible benchmark that doesn't fabricate data or
report invalid metrics.
"""