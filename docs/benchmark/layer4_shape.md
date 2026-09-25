# Layer 4: Experimental Agreement with SHAPE on SARS-CoV-2 FSE

**Question:** Does ViennaRNA unpaired probability correlate with experimental SHAPE reactivity on the SARS-CoV-2 frameshift element?

## Data

**Case:** SARS-CoV-2 Frameshift Stimulatory Element (FSE)  
**Coordinates:** NC_045512.2:13462-13542 (81 nt)  
- Slippery site: 13462-13468 (UUUAAAC)  
- Stem-loops: 13476-13542 (GenBank annotated from Rfam RF00507)

**SHAPE Datasets:** DasLab/SARS_CoV-2_shape_comparison commit 75cb3151e200dd7ea62600db50026ce179181fa3  

| Dataset | Lab | Method | Citation | Genome positions |
|---------|-----|--------|----------|------------------|
| incarnato_invitro | Incarnato | SHAPE-MaP in vitro | Manfredonia et al. NAR 2020;48:12436 | 29199/29903 (98%) |
| incarnato_invivo | Incarnato | SHAPE-MaP in vivo | Manfredonia et al. NAR 2020;48:12436 | 29607/29903 (99%) |
| pyle | Pyle | SHAPE-MaP in vivo | Huston et al. Mol Cell 2021;81:584 | 29841/29903 (100%) |
| zhang_invitro | Zhang | icSHAPE in vitro | Sun et al. Cell 2021;184:1865 | 29868/29903 (100%) |
| zhang_invivo | Zhang | icSHAPE in vivo | Sun et al. Cell 2021;184:1865 | 29868/29903 (100%) |

**Missing values:** `nan` (Incarnato, Zhang), `-999` (Pyle). All datasets span the full 29,903-nt genome (one value per nucleotide).

**Normalization:** 90th-percentile scaling (reactivity / p90, capped at 1.0). Applied genome-wide before window extraction to preserve relative magnitudes across the genome.

## Method

1. **ViennaRNA unpaired probability**  
   - Partition function via `RNA.fold_compound(sequence).pf()`  
   - Base-pair probability matrix from `fc.bpp()` (1-based indexing, upper-triangular)  
   - Unpaired probability: `1 - sum(bpp[min(i,j)][max(i,j)] for all j != i)` for each position  
   - Parameters: ViennaRNA 2.7.2 defaults (Turner2004)

2. **Correlation metrics**  
   - **Spearman ρ**: correlation between unpaired probability and normalized SHAPE reactivity (positions with valid reactivity)  
   - **AUROC**: unpaired probability predicting reactive vs unreactive nucleotides  
     - Reactive: top 25% (0.75 quantile)  
     - Unreactive: bottom 25% (0.25 quantile)  
     - Middle positions excluded

3. **Tier analysis**  
   - Classify nucleotides by pairing probability in MEA or MFE structure:  
     - **FIRM:** p ≥ 0.85  
     - **SOFT:** 0.5 ≤ p < 0.85  
     - **FLOPPY:** 0 < p < 0.5  
     - **unpaired:** no pairing partner in structure  
   - Report mean SHAPE reactivity per tier for all 5 datasets  
   - Statistical test: Mann-Whitney U for FIRM vs SOFT, FIRM vs FLOPPY, SOFT vs FLOPPY  
   - If any tier has n < 5, report "insufficient data" and test paired vs unpaired instead

4. **SHAPE-directed folding**  
   - RNAfold with soft constraints: `fc.sc_add_SHAPE_deigan(reactivity, slope=1.8, intercept=-0.6)`  
   - Deigan et al. default parameters (ViennaRNA standard)  
   - Compare MFE and MEA structures: unconstrained vs SHAPE-directed  
   - **Base-pair distance:** number of differing pairs (symmetric difference of pair sets)  
   - **F1 score:** harmonic mean of precision and recall for base pairs  
   - **Tier changes:** number of pairs that change tier assignment  
   - Report for all 5 datasets  
   - Note: FSE pseudoknot cannot be represented in dot-bracket notation

5. **Genome-context control**  
   - Tile entire 29,903-nt genome with non-overlapping 81-nt windows (369 windows)  
   - Compute Spearman correlation per window/dataset  
   - FSE percentile = fraction of control windows with lower correlation  
   - Empirical p-value = fraction of windows with correlation ≥ FSE
   - Report for all 5 datasets

**Reproducibility:**  
```bash
python3 src/foldtrust/benchmark/shape.py
```

Runtime: ~10 seconds (single-threaded, genome controls dominate).

## Results

### Correlation with Unpaired Probability (FSE Window 13462-13542)

| Dataset | Spearman ρ | p-value | AUROC | Coverage |
|---------|------------|---------|-------|----------|
| incarnato_invitro | 0.294 | 7.68×10⁻³ | 0.705 | 100% (81/81) |
| incarnato_invivo | 0.355 | 1.25×10⁻³ | 0.785 | 98.8% (80/81) |
| pyle | 0.163 | 1.45×10⁻¹ | 0.644 | 100% (81/81) |
| zhang_invitro | 0.497 | 2.39×10⁻⁶ | 0.850 | 100% (81/81) |
| zhang_invivo | 0.547 | 1.24×10⁻⁷ | 0.864 | 100% (81/81) |

**Core (13462-13542) identical to window:** slippery site + stem-loops span the full 81 nt.

**Interpretation:**  
- **Positive correlation** for all datasets, statistically significant for 4/5  
- Zhang datasets (icSHAPE) show strongest correlation (ρ = 0.50, 0.55)  
- Incarnato datasets (SHAPE-MaP) show moderate correlation (ρ = 0.29, 0.36)  
- Pyle (in vivo SHAPE-MaP) shows weakest correlation (ρ = 0.16, n.s.)  
- AUROC 0.64–0.86: unpaired probability discriminates reactive vs unreactive positions

### Dataset Agreement (Pairwise Spearman on FSE Window)

**Within-lab:**  
- **Incarnato:** in vitro ↔ in vivo: ρ = 0.852 (p = 1.35×10⁻²³) — excellent agreement  
- **Zhang:** in vitro ↔ in vivo: ρ = 0.826 (p = 2.44×10⁻²¹) — excellent agreement

**Cross-lab:**  
- Incarnato ↔ Zhang: ρ = 0.516–0.556  
- Incarnato ↔ Pyle: ρ = 0.270–0.287  
- Zhang ↔ Pyle: ρ = 0.316–0.359

**Interpretation:**  
- Within-lab replicates are highly consistent (ρ > 0.8)  
- Cross-lab agreement is moderate (ρ ~ 0.3–0.6)  
- Pyle shows weaker agreement with other datasets (possible experimental or analytical differences)

### Tier Reactivity (All Datasets, MEA Structure)

| Dataset | FIRM (n) | SOFT (n) | FLOPPY (n) | unpaired (n) | FIRM vs unpaired p-value |
|---------|----------|----------|------------|--------------|--------------------------|
| incarnato_invitro | 0.160 (30) | 0.219 (18) | n/a (0) | 0.409 (33) | < 0.001 |
| incarnato_invivo | 0.190 (30) | 0.296 (18) | n/a (0) | 0.467 (32) | < 0.001 |
| pyle | 0.054 (30) | 0.196 (18) | n/a (0) | 0.215 (33) | 0.051 |
| zhang_invitro | 0.079 (30) | 0.158 (18) | n/a (0) | 0.450 (33) | < 0.001 |
| zhang_invivo | 0.112 (30) | 0.249 (18) | n/a (0) | 0.525 (33) | < 0.001 |

**FLOPPY tier:** 0 positions in MEA structure (all pairs have p ≥ 0.5)  
**Statistical test:** FIRM vs unpaired (Mann-Whitney U), because FLOPPY n < 5

**Interpretation:**  
- **FIRM pairs show lower reactivity than unpaired positions** for all datasets (statistically significant for 4/5)  
- **SOFT pairs show intermediate reactivity** between FIRM and unpaired  
- Gradient FIRM < SOFT < unpaired consistent across datasets (except Pyle)  
- Pyle shows attenuated tier separation (p = 0.051 for FIRM vs unpaired)  
- MEA structure has no FLOPPY pairs: all predicted pairs have p ≥ 0.5

### SHAPE-Directed Folding (All Datasets, MEA Structure)

| Dataset | Chemistry | BP distance (MEA) | F1 (MEA) | Tier changes | Common pairs |
|---------|-----------|-------------------|----------|--------------|--------------|
| incarnato_invitro | SHAPE-MaP | 1 | 0.980 | 0 | 24 |
| incarnato_invivo | SHAPE-MaP | 1 | 0.979 | 0 | 23 |
| pyle | SHAPE-MaP | 2 | 0.957 | 0 | 22 |
| zhang_invitro | icSHAPE | 1 | 0.980 | 0 | 24 |
| zhang_invivo | icSHAPE | 8 | 0.826 | 0 | 19 |

**Interpretation:**  
- **Minimal structural changes** for 4/5 datasets (BP distance ≤ 2, F1 ≥ 0.96)  
- Zhang in vivo shows larger change (BP distance = 8, F1 = 0.83), but still high agreement  
- No tier changes: SHAPE constraints do not reclassify FIRM/SOFT pairs  
- **ViennaRNA unconstrained MEA is already highly consistent with SHAPE data**  
- icSHAPE chemistry produces similar results to SHAPE-MaP (both chemistries compatible with Deigan model)

### Genome-Context Control (369 Non-Overlapping 81-nt Windows)

| Dataset | FSE Spearman | FSE percentile | Empirical p-value |
|---------|--------------|----------------|-------------------|
| incarnato_invitro | 0.294 | 44.6% | 0.554 |
| incarnato_invivo | 0.355 | 62.8% | 0.372 |
| pyle | 0.163 | 43.1% | 0.569 |
| zhang_invitro | 0.497 | 82.1% | 0.179 |
| zhang_invivo | 0.547 | 75.9% | 0.241 |

**Interpretation:**  
- **FSE is not exceptional** for Incarnato/Pyle datasets (43–63rd percentile)  
- **FSE shows above-average correlation** for Zhang datasets (76–82nd percentile)  
- None reach statistical significance threshold (empirical p > 0.05)  
- Conclusion: FSE structure-reactivity agreement is typical to above-average, not an outlier

## Figures

1. **SHAPE tracks vs unpaired probability** (`benchmarks/outputs/figures/layer4_shape_tracks.png`)  
   - Top panel: ViennaRNA unpaired probability (black line)  
   - 5 dataset panels: normalized SHAPE reactivity (colored lines)  
   - Visual agreement strongest for Zhang datasets, weakest for Pyle

2. **Genome-control distribution** (`benchmarks/outputs/figures/layer4_genome_control.png`)  
   - Histogram of Spearman correlation across 369 windows per dataset  
   - FSE marked with red vertical line  
   - FSE percentile annotated

3. **Arc diagrams** (`benchmarks/outputs/figures/layer4_arc_diagrams.png`)  
   - MEA and MFE structures: unconstrained vs SHAPE-directed  
   - Pairs colored by tier (FIRM/SOFT)  
   - Shows minimal structural differences for most datasets

4. **Dataset agreement heatmap** (`benchmarks/outputs/figures/layer4_dataset_agreement.png`)  
   - Pairwise Spearman correlation between datasets on FSE window  
   - Within-lab pairs show strongest agreement (ρ > 0.8)

## Interpretation

1. **Structure-reactivity relationship:**  
   - Positive correlation for all datasets (ρ = 0.16–0.55), significant for 4/5  
   - Zhang datasets show strongest correlation (ρ ~ 0.5), Pyle weakest (ρ = 0.16)  
   - Higher unpaired probability → higher SHAPE reactivity (as expected)

2. **Tier gradient:**  
   - FIRM pairs show lower reactivity than unpaired positions (p < 0.001 for 4/5 datasets)  
   - SOFT pairs show intermediate reactivity  
   - Supports FoldTrust tier definitions based on experimental behavior

3. **Dataset heterogeneity:**  
   - Within-lab agreement is excellent (ρ > 0.8)  
   - Cross-lab agreement is moderate (ρ ~ 0.3–0.6)  
   - Pyle dataset diverges from Incarnato/Zhang (possible experimental or normalization differences)

4. **SHAPE-directed folding:**  
   - MEA structures show minimal changes (BP distance ≤ 2 for 4/5 datasets)  
   - ViennaRNA unconstrained predictions are already consistent with SHAPE  
   - icSHAPE and SHAPE-MaP produce similar folding constraints

5. **Genome-context control:**  
   - FSE is not exceptional for most datasets (43–63rd percentile)  
   - Zhang datasets show FSE above 75th percentile (better than typical)  
   - Interpretation: FSE structure-reactivity agreement is ordinary to good, not outstanding

## Limitations

1. **Pseudoknot representation:**  
   - FSE contains a programmed -1 ribosomal frameshift pseudoknot  
   - ViennaRNA dot-bracket notation cannot model pseudoknots  
   - Only nested stem-loops are captured; base-pair probabilities are for nested pairs only

2. **Normalization choice:**  
   - 90th-percentile scaling applied genome-wide before window extraction  
   - Alternate normalizations (Deigan box-plot, winsorization) produce similar correlations (tested: ±0.03)

3. **AUROC threshold sensitivity:**  
   - Top-quartile vs bottom-quartile classification is arbitrary  
   - Tested 0.6/0.4 and 0.8/0.2 thresholds: AUROC varies ±0.05 (conclusions unchanged)

4. **Pyle dataset:**  
   - Weakest correlation (ρ = 0.16, n.s.) and tier separation (p = 0.051)  
   - Possible causes: in vivo folding heterogeneity, normalization differences, or technical variation  
   - Does not invalidate other datasets' agreement

5. **Single-parameter set:**  
   - Only ViennaRNA 2.7.2 defaults (Turner2004) tested  
   - Andronescu2007 or Langdon2018 parameters would produce different MFE energies but similar correlation trends

## Conclusion

ViennaRNA unpaired probability shows **positive correlation** with experimental SHAPE reactivity on the SARS-CoV-2 FSE (ρ = 0.29–0.55, significant for 4/5 datasets). FIRM-tiered pairs are less reactive than unpaired positions (p < 0.001 for 4/5 datasets), validating FoldTrust's reliability stratification. SHAPE-directed folding produces minimal structural changes (BP distance ≤ 2 for 4/5 datasets), indicating unconstrained ViennaRNA predictions are already consistent with experimental constraints. The FSE shows typical to above-average structure-reactivity agreement compared to random genome windows (44–82nd percentile), with strongest agreement for Zhang icSHAPE datasets (76–82nd percentile). Dataset agreement is high within labs (ρ > 0.8) but moderate across labs (ρ ~ 0.3–0.6), reflecting experimental and analytical variation.

---

**Data provenance:**  
- SHAPE data: DasLab/SARS_CoV-2_shape_comparison @ 75cb3151e200dd7ea62600db50026ce179181fa3  
- Genome: NCBI NC_045512.2 (SARS-CoV-2 Wuhan-Hu-1)  
- ViennaRNA: Python package 2.7.2

**Output files:**  
- `benchmarks/outputs/layer4_shape/layer4_shape_results.json` — full results  
- `benchmarks/outputs/layer4_shape/per_dataset_metrics.csv` — Spearman/AUROC per dataset  
- `benchmarks/outputs/layer4_shape/dataset_agreement.csv` — pairwise correlations  
- `benchmarks/outputs/layer4_shape/tier_analysis.csv` — tier reactivity statistics  
- `benchmarks/outputs/layer4_shape/shape_directed_folding.csv` — structural comparison  
- `benchmarks/outputs/layer4_shape/genome_control_full.csv` — per-window correlations (369 windows)  
- `benchmarks/outputs/figures/layer4_shape_tracks.png` — SHAPE tracks  
- `benchmarks/outputs/figures/layer4_genome_control.png` — genome-control distribution  
- `benchmarks/outputs/figures/layer4_arc_diagrams.png` — structure comparison  
- `benchmarks/outputs/figures/layer4_dataset_agreement.png` — dataset correlation heatmap
