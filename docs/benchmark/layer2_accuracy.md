# Layer 2: Structure Accuracy

## Question

How accurate are ViennaRNA MFE, MEA, and centroid structures compared to known reference structures?

## Data

Three reference datasets, filtered to sequences ≤500 nucleotides:

- **ArchiveII** (Mathews lab): 3,854 structures (filtered from 3,975 total)
  - Source: Comparative analysis and crystallography
  - Citation: Sloma MF, Mathews DH. RNA 2016;22:1808-1818. DOI [10.1261/rna.053694.115](https://doi.org/10.1261/rna.053694.115)
  - Families: 16S rRNA, 23S rRNA, 5S rRNA, RNase P, group I/II introns, SRP, tRNA, telomerase RNA, tmRNA

- **bpRNA TS0** (test set 0): 1,305 structures
  - Source: Comparative analysis (CRW, tmRNA website, Rfam)
  - Citation: MXfold2 Zenodo release [10.5281/zenodo.4430150](https://doi.org/10.5281/zenodo.4430150)
  - Length range: 22–499 nt (mean 136.1 nt)

- **Rfam seed** (release 15.1): 468 structures
  - Source: Rfam consensus structure projected onto seed sequences
  - Citation: Ontiveros-Palacios N, et al. Nucleic Acids Res 2025;53(D1):D258-D267. DOI [10.1093/nar/gkae1023](https://doi.org/10.1093/nar/gkae1023)
  - 19 families: tRNA, 5S rRNA, SRP, RNase P, spliceosomal RNAs (U1/U2/U4/U6), riboswitches (FMN, SAM, purine, TPP, cobalamin, glycine, magnesium, lysine, glmS), tmRNA

**Sampling**: Random seed 42, 200 structures per dataset (re-seeded before each draw to ensure reproducibility). Total sample: 600 structures.

**Pseudoknot handling**: Reference structures with crossing pairs (pseudoknots) are reduced to nested-only structure using greedy selection: pairs are sorted by 5' position (ties broken by longer span) and added if they do not cross any already-selected pair.

**Known limitations**: 
- Nine duplicate-sequence pairs exist in the 600-structure sample (one pair spans ArchiveII/Rfam; not deduplicated)
- ArchiveII includes both canonical and non-canonical pairs; ViennaRNA predicts only canonical/wobble pairs

## Method

For each structure:
1. Fold sequence with ViennaRNA 2.7.2 (Turner2004 parameters, 37°C)
2. Predict three structures:
   - **MFE**: Minimum free energy structure (`RNA.fold_compound.mfe()`)
   - **MEA**: Maximum expected accuracy (γ=1, `RNA.fold_compound.MEA()`)
   - **Centroid**: Centroid structure (`RNA.fold_compound.centroid()`)
3. Compare predicted pairs to reference pairs (after PK removal)
4. Compute per-structure sensitivity, PPV, F1, MCC
   - **MCC** (Matthews Correlation Coefficient): MCC = sqrt(sensitivity × PPV), the standard approximation for RNA secondary structure (Gorodkin, Stricklin & Stormo 2001, *Nucleic Acids Res* 29:2135-2144, [doi:10.1093/nar/29.10.2135](https://doi.org/10.1093/nar/29.10.2135))
5. Report mean of per-structure F1 (primary metric)

**Matching conventions**:
- **Exact matching** (primary): Predicted pair (i,j) is correct if (i,j) is in the reference
- **Slip-tolerant matching**: Predicted pair (i,j) matches if any of (i,j), (i±1,j), (i,j±1) is in the reference. A reference pair is recovered if any of (i,j), (i±1,j), (i,j±1) is predicted. By construction, slip F1 ≥ exact F1.

**Averaging**:
- Primary metric: mean of per-structure F1
- Pooled metric: F1 computed from sum of all TP, FP, FN across structures
- Bootstrap 95% confidence intervals (1,000 resamples)

**Command**:
```bash
python scripts/run_layers_1_2_3.py
```

## Results

### Overall Metrics (Mean of Per-Structure)

| Method | Sensitivity | PPV | F1 | MCC |
|--------|-------------|-----|----|----|
| MFE | 0.639 | 0.499 | 0.548 | 0.557 |
| MEA | 0.641 | 0.521 | 0.563 | 0.571 |
| CENTROID | 0.622 | 0.551 | 0.570 | 0.577 |

### Per-Dataset F1 (Mean of Per-Structure)

| Dataset | N | MFE | MEA | Centroid |
|---------|---|----|-----|---------
| ArchiveII | 200 | 0.562 | 0.582 | 0.591 |
| Rfam_seed | 200 | 0.569 | 0.582 | 0.586 |
| bpRNA_TS0 | 200 | 0.514 | 0.525 | 0.532 |

**Interpretation**:

- **Centroid is the most accurate** single structure (F1 0.570), followed by MEA (0.563) and MFE (0.548)
- **bpRNA TS0 is hardest** (F1 ~0.52), likely due to higher structural diversity and fewer constraints from covariation
- **ArchiveII and Rfam are similar** (F1 ~0.58), benefiting from strong comparative-analysis signals
- All methods show moderate to good correlation (MCC 0.557–0.577)

## Limitations

1. **Reference bias**: Comparative structures may already reflect thermodynamic preferences that ViennaRNA captures
2. **Pseudoknot removal**: Some functional pairs are discarded from references
3. **Non-canonical pairs**: ArchiveII includes non-canonical pairs that ViennaRNA cannot predict
4. **Sample size**: 200 per dataset is sufficient for means but confidence intervals overlap
5. **Length filter**: Excludes long RNAs (>500 nt) where ViennaRNA accuracy may differ
6. **Turner2004 only**: Other parameter sets (Andronescu2007, Langdon2018) not tested in this layer
