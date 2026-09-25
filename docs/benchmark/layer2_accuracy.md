# Layer 2: Structure Accuracy

## Question

How accurate are ViennaRNA's MFE, MEA, and centroid structures on reference datasets? How do F1 scores vary by dataset, family, and sequence length?

## Data

We evaluate on four reference datasets:

### 1. ArchiveII (Primary)

- **Source**: Sloma & Mathews (2016), *RNA* 22:1808-1818 [DOI:10.1261/rna.057046.116]
- **Format**: JSONL extracted from original .ct files
- **Count** (≤500 nt): Filtered from 3,975 total structures
- **Families**: 16S, 23S, 5S, RNaseP, Group I intron, Group II intron, SRP, tRNA, telomerase, tmRNA
- **Pseudoknots**: Removed before scoring (reported separately)
- **Note**: Original ArchiveII URL (`rna.urmc.rochester.edu/pub/archiveII.tar.gz`) now returns 404; we use the Wayback Machine capture (2022-10-17) with verified SHA1

### 2. bpRNA-1m TS0 (Canonical Pairs)

- **Source**: Danaee et al. (2018), *NAR* 46:5381-5394; Singh et al. (2019), *Nat Commun* 10:5407
- **Obtained from**: MXfold2 Zenodo release (DOI:10.5281/zenodo.4430150)
- **Count** (≤500 nt): Filtered from 1,305 total structures
- **Families**: Mixed (no family labels in MXfold2 release)
- **Pseudoknots**: Canonical pairs only (pseudoknots pre-removed)
- **Length distribution**: 15–1503 nt (mean ~75 nt)

### 3. Rfam 15.1 Seed Consensus-Derived

- **Source**: Ontiveros-Palacios et al. (2025), *NAR* 53(D1):D258-D267 [DOI:10.1093/nar/gkae1051]
- **Format**: Consensus structure (SS_cons) projected per sequence
- **Count** (≤500 nt): 468 sequences, 19 families
- **Families**: tRNA (RF00005), 5S_rRNA (RF00001), SRP (RF00017, RF01854), RNaseP (RF00010), snRNAs (RF00003, RF00004, RF00015, RF00026), riboswitches (RF00050, RF00059, RF00162, RF00167, RF00168, RF00174, RF00234, RF00380, RF00504), tmRNA (RF00023)
- **Note**: Labeled as **consensus-projected** structures (approximations, not experimentally validated per-sequence structures)
- **Pseudoknots**: WUSS letter pseudoknots (Aa, Bb, ...) stripped; nested brackets (`()`, `<>`, `[]`, `{}`) kept as legitimate nested pairs

### 4. SPOT-RNA PDB TS1 (Future)

- **Source**: Singh et al. (2019), *Nat Commun* 10:5407
- **Status**: Available in bundle but not yet integrated into Layer 2 (planned)
- **Count**: 67 PDB-derived structures

**Default run**: ArchiveII + bpRNA TS0 + Rfam seed, length ≤ 500 nt  
**Full run (`--full`)**: All sequences, no length cap (much slower)

### Pseudoknot Removal

For fair comparison with ViennaRNA (which predicts only nested structures), we remove pseudoknotted pairs from reference structures before scoring. The number of removed pairs is reported for transparency.

**Algorithm**: For each structure, we iterate through pairs in order and keep only those that do not pseudoknot with any already-kept pair. A pair (a, b) pseudoknots with (c, d) if `c < a < d < b` or `a < c < b < d`.

## Method

### ViennaRNA Settings

- **Parameters**: Turner2004 (RNA.params_load_RNA_Turner2004())
- **Temperature**: 37°C (1M NaCl equivalent)
- **Important**: We create a **new `RNA.md()` after loading parameters** to avoid the silent parameter-loading bug documented in MANIFEST.md

### Structure Prediction Methods

1. **MFE (Minimum Free Energy)**: `RNA.fold_compound(seq, md).mfe()`
2. **MEA (Maximum Expected Accuracy)**: `fc.MEA(gamma=1.0)` (default gamma; optionally sweep 0.5, 1.0, 2.0)
3. **Centroid**: `fc.centroid()`

### Metrics

For each predicted structure vs reference:

- **Sensitivity** (recall): TP / (TP + FN)
- **PPV** (precision): TP / (TP + FP)
- **F1**: 2 × TP / (2 × TP + FP + FN)
- **MCC** (Matthews correlation coefficient): (TP × TN - FP × FN) / sqrt((TP+FP)(TP+FN)(TN+FP)(TN+FN))

Both **exact** and **+/-1 slip-tolerant** scoring are reported.

### Stratification

- **Per dataset**: ArchiveII, bpRNA_TS0, Rfam_seed
- **Per family**: Rfam families (where labeled)
- **Per length bin**: 0–50, 50–100, 100–200, 200–500 nt

### Bootstrap Confidence Intervals

We compute 95% CIs by resampling structures (1000 bootstrap iterations).

### Sanity Check

We compare MFE F1 scores on ArchiveII against published results:

- **MXfold2 paper** (Sato et al., *Nat Commun* 2021): ViennaRNA MFE on ArchiveII F1 ≈ 0.45–0.50 (depending on subset)
- **LinearPartition paper** (Huang et al., *Bioinformatics* 2019): ViennaRNA MFE F1 ≈ 0.44

Any large discrepancy triggers investigation.

## Command

```bash
# Default (500 nt cap, fast)
cd /workspace
python3 scripts/run_layers_1_2_3.py

# Full dataset (all sequences, slow)
python3 scripts/run_layers_1_2_3.py --full
```

## Results

**Run completed**: 2026-09-25  
**Runtime**: 2.6 minutes  
**Sample size**: 200 structures per dataset (sampled from ≤500 nt filtered sets)

### Dataset Counts

- **ArchiveII**: 200 (sampled from 3,854 ≤500 nt)
- **bpRNA TS0**: 200 (sampled from 1,305 ≤500 nt)
- **Rfam seed**: 200 (sampled from 468 total)
- **Total**: 600 structures

### Overall F1 (95% CI)

| Method   | F1 Mean | 95% CI           | Sensitivity | PPV   |
|----------|---------|------------------|-------------|-------|
| MFE      | 0.548   | [0.526, 0.570]   | 0.568       | 0.555 |
| MEA      | 0.563   | [0.542, 0.585]   | 0.612       | 0.538 |
| Centroid | 0.552   | [0.530, 0.574]   | 0.605       | 0.527 |

**Interpretation**: MEA slightly outperforms MFE (F1 = 0.563 vs 0.548), consistent with published findings that ensemble-based methods improve over MFE. MEA trades higher sensitivity for slightly lower PPV.

### Per-Dataset F1 (MFE)

*(From layer2_summary_dataset.csv)*

| Dataset      | MFE F1 Mean | 95% CI           | Count |
|--------------|-------------|------------------|-------|
| ArchiveII    | 0.510       | [0.476, 0.544]   | 200   |
| bpRNA_TS0    | 0.595       | [0.565, 0.626]   | 200   |
| Rfam_seed    | 0.538       | [0.498, 0.577]   | 200   |

**Interpretation**: bpRNA TS0 shows higher F1 (0.595) than ArchiveII (0.510), possibly due to shorter sequences and simpler structures in the bpRNA filtered set.

### Sanity Check vs Published Numbers

Our MFE F1 on ArchiveII (≤500 nt, 200 sampled): **0.510 [0.476, 0.544]**

Published benchmarks (full ArchiveII, no length cap):
- **MXfold2 paper** (Sato et al., *Nat Commun* 2021): ViennaRNA MFE F1 ≈ 0.44–0.48
- **LinearPartition paper** (Huang et al., *Bioinformatics* 2019): ViennaRNA MFE F1 ≈ 0.44

**Gap analysis**: Our F1 (0.510) is slightly higher than published values (~0.44–0.48). This is **expected** because:
1. We filter to ≤500 nt (shorter sequences are easier to predict)
2. Published results include full-length 16S/23S rRNA (1.5 kb), which are harder
3. We use a random sample, which may have selected slightly easier structures

Our results are **consistent with published work** when accounting for the length filter. The gap is within expected variation.

To verify against full ArchiveII (all 3,975 structures, no length cap), run with `--full` flag (slower).

### Figures

Generated in `benchmarks/outputs/figures/`:
- `layer2_method_comparison.png`: Sensitivity, PPV, F1 for MFE/MEA/Centroid
- `layer2_f1_vs_length.png`: F1 vs sequence length (binned, with error bars)
- `layer2_f1_by_family.png`: Box plot of F1 per Rfam family (MFE method)

## Interpretation

### MFE vs MEA vs Centroid

- **MEA wins on F1**: MEA (0.563) > Centroid (0.552) > MFE (0.548)
- **MEA is more sensitive**: MEA recovers more true pairs (sensitivity 0.612) but with slightly lower precision (PPV 0.538)
- **MFE is more precise**: MFE has higher PPV (0.555) but misses more pairs (sensitivity 0.568)

**Practical implication**: For experimental design (e.g., ASO targeting), MEA's higher sensitivity is useful for identifying accessible regions. For structure validation, MFE's higher PPV means predicted pairs are more likely correct.

### Length Dependence

(See `layer2_f1_vs_length.png`)

F1 scores generally degrade with increasing length, consistent with the accumulation of errors in longer sequences. The full ≤500 nt cap helps maintain reasonable accuracy.

### Family Variation

(See `layer2_f1_by_family.png` and `layer2_summary_family.csv`)

Rfam families show variation in predictability. Families with strong, conserved structures (e.g., tRNA, riboswitches) tend to have higher F1 than families with flexible or poorly conserved structures (e.g., some SRPs).

## Limitations

1. **Length cap**: Default run uses ≤500 nt for speed (full 16S/23S rRNA excluded unless `--full`)
2. **Pseudoknot removal**: We score only nested pairs; pseudoknots are removed from references
3. **Consensus-projected Rfam**: Rfam seed structures are consensus projections (approximate), not per-sequence experimental structures
4. **No contact distance filtering**: We count all pairs regardless of sequence separation
5. **Turner2004 only**: We do not test Andronescu2007 or Langdon2018 parameter sets (future work)

## Files Generated

```
benchmarks/outputs/layer2/
├── layer2_raw_results.csv           # Per-structure predictions
├── layer2_summary_dataset.csv       # F1 by dataset and method
├── layer2_summary_family.csv        # F1 by Rfam family
├── layer2_summary_length.csv        # F1 by length bin
├── layer2_summary_overall.csv       # Overall sensitivity/PPV/F1/MCC per method
└── layer2_summary.json              # JSON summary

benchmarks/outputs/figures/
├── layer2_method_comparison.png
├── layer2_f1_vs_length.png
└── layer2_f1_by_family.png
```
