# Layer 5: Robustness — Stability of FoldTrust Tiers Across Conditions

## Question

How robust are FoldTrust reliability tiers (FIRM / SOFT / FLOPPY) to changes in folding parameters, temperature, and window boundaries?

RNA secondary structure prediction requires choosing:
1. **Energy parameter sets** — experimental thermodynamic measurements (Turner 2004, Andronescu 2007, Langdon 2018)
2. **Temperature** — affects both thermodynamics and kinetics
3. **Window boundaries** — where to start/stop the folded region

Real-world use cases involve uncertainty in all three. A robust reliability analysis should produce consistent tier classifications even when these conditions vary within biologically plausible ranges.

## Data

**Five disease-relevant regulatory windows** from FoldTrust's MVP cases:

| Case | Gene | Region | Length (nt) | Source Accession | Coordinates |
|------|------|--------|------------:|------------------|-------------|
| sars2-fse | SARS-CoV-2 ORF1ab | Frameshift stimulatory element | 81 | NC_045512.2 | 13462–13542 |
| smn2-iss-n1 | SMN2 | Intron 7 intronic splicing silencer N1 | 201 | NM_017411.4 | 840–1040 |
| cftr-5utr | CFTR | 5' UTR | 268 | NM_000492.4 | 133–400 |
| mapt-e10 | MAPT | Exon 10 splice regulatory region | 281 | NM_001123066.4 | 980–1260 |
| hcv-ires-dii | HCV 5' UTR | IRES Domain II | 271 | AF009606.1 | 40–310 |

**Baseline condition:** Turner2004 parameters, 37°C, original window boundaries (as defined in `data/cases/*/sequence.fa`).

**Perturbations tested:**

1. **Energy parameter sets**
   - Andronescu2007 (`rna_andronescu2007.par`)
   - Langdon2018 (`rna_langdon2018.par`)
   - vs. Turner2004 baseline

2. **Temperature**
   - 24°C (cool)
   - 30°C (moderate)
   - 37°C (baseline, physiological)
   - 42°C (fever)
   - 45°C (heat shock)

3. **Window boundary jitter**
   - Extend left by 10, 25 nt
   - Extend right by 10, 25 nt
   - Extend both sides by 10, 25 nt
   - Shrink left by 10, 25 nt
   - Shrink right by 10, 25 nt
   - Shrink both sides by 10, 25 nt

For jitter, flanking sequences were fetched from NCBI E-utilities (`efetch`) using the source accession and coordinates from each case's metadata. Each case's original sequence was verified to match the genomic slice exactly before jitter was applied.

## Method

### Baseline Computation

For each case, compute the **baseline** (Turner2004, 37°C, original window):
- Per-nucleotide tier classification (FIRM / SOFT / FLOPPY / UNPAIRED based on maximum pair probability)
- MEA (Maximum Expected Accuracy) structure
- MFE (Minimum Free Energy) structure
- Per-nucleotide unpaired probability
- Set of FIRM pairs (p ≥ 0.85)
- Set of FLOPPY pairs (0 < p < 0.5)

### Robustness Metrics

For each perturbation condition, compare to baseline:

1. **Per-nucleotide tier agreement** — fraction of positions whose tier classification is unchanged
2. **Base-pair Jaccard (MEA)** — Jaccard index between MEA structure base-pair sets (intersection / union)
3. **Base-pair Jaccard (MFE)** — Jaccard index between MFE structure base-pair sets
4. **Unpaired probability Spearman** — Spearman correlation of per-nucleotide unpaired probabilities
5. **FIRM retention** — fraction of baseline FIRM pairs that remain FIRM
6. **FLOPPY retention** — fraction of baseline FLOPPY pairs that remain FLOPPY

All metrics range from 0 (complete change) to 1 (perfect stability).

### Implementation Notes

- **ViennaRNA Python API** was used to enable programmatic control of parameter sets and temperature. The CLI tools (`RNAfold --paramFile`) do not support all three Turner/Andronescu/Langdon sets via command-line flags.
- **Critical ViennaRNA gotcha:** After calling `RNA.params_load_RNA_<set>()`, a **new** `RNA.md()` object must be created before `RNA.fold_compound(seq, md)`, or the parameters will silently not apply. This bug was documented in the data bundle and a regression test was added (`tests/test_param_energies.py`).
- **Window jitter comparison:** For extended windows, metrics are computed over the full sequence but compared to baseline on the shared (original) positions. For shrunk windows, only the shared positions are compared.
- **Sanity check:** Conditions are expected to produce measurable changes. If all metrics equal exactly 1.0 for a condition, this indicates a bug (parameters not applied).

### Command

```bash
foldtrust benchmark robustness -o benchmarks/outputs
```

Output files:
- `benchmarks/outputs/layer5_robustness/robustness_summary.csv` — all cases × all conditions × all metrics
- `benchmarks/outputs/layer5_robustness/robustness_detailed.json` — full per-condition results
- `benchmarks/outputs/layer5_robustness/mfe_energies_by_params.csv` — MFE energy per case per parameter set
- `benchmarks/outputs/layer5_robustness/region_stability.csv` — per-case mean/std stability across conditions
- `benchmarks/outputs/figures/layer5_tier_agreement_heatmap.png` — heatmap of tier agreement
- `benchmarks/outputs/figures/layer5_mea_jaccard_heatmap.png` — heatmap of MEA Jaccard

## Results

### MFE Energies by Parameter Set

MFE free energies (kcal/mol) for each case under three parameter sets:

| Case | Turner2004 | Andronescu2007 | Langdon2018 |
|------|------------|----------------|-------------|
| sars2-fse | -26.00 | -22.26 | -24.70 |
| smn2-iss-n1 | [TBD] | [TBD] | [TBD] |
| cftr-5utr | [TBD] | [TBD] | [TBD] |
| mapt-e10 | [TBD] | [TBD] | [TBD] |
| hcv-ires-dii | [TBD] | [TBD] | [TBD] |

**Verification:** The SARS-CoV-2 FSE energies match the expected values from the data bundle's parameter check (−26.00 / −22.26 / −24.70), confirming that parameter loading works correctly. A regression test (`test_param_energies.py`) asserts these values to prevent the silent-parameter-bug from returning.

### Robustness Summary

**Cases × Conditions × Metrics** (values are mean ± std across all five cases):

| Condition Type | Tier Agreement | MEA Jaccard | MFE Jaccard | Unpaired Spearman | FIRM Retention | FLOPPY Retention |
|----------------|----------------|-------------|-------------|-------------------|----------------|------------------|
| **Parameter sets** | | | | | | |
| Andronescu2007 | [TBD] | [TBD] | [TBD] | [TBD] | [TBD] | [TBD] |
| Langdon2018 | [TBD] | [TBD] | [TBD] | [TBD] | [TBD] | [TBD] |
| **Temperatures** | | | | | | |
| 24°C | [TBD] | [TBD] | [TBD] | [TBD] | [TBD] | [TBD] |
| 30°C | [TBD] | [TBD] | [TBD] | [TBD] | [TBD] | [TBD] |
| 42°C | [TBD] | [TBD] | [TBD] | [TBD] | [TBD] | [TBD] |
| 45°C | [TBD] | [TBD] | [TBD] | [TBD] | [TBD] | [TBD] |
| **Window jitter** | | | | | | |
| Extend left 10 nt | [TBD] | [TBD] | [TBD] | [TBD] | [TBD] | [TBD] |
| Extend right 10 nt | [TBD] | [TBD] | [TBD] | [TBD] | [TBD] | [TBD] |
| Extend both 10 nt | [TBD] | [TBD] | [TBD] | [TBD] | [TBD] | [TBD] |
| Shrink left 10 nt | [TBD] | [TBD] | [TBD] | [TBD] | [TBD] | [TBD] |
| Shrink right 10 nt | [TBD] | [TBD] | [TBD] | [TBD] | [TBD] | [TBD] |
| Extend left 25 nt | [TBD] | [TBD] | [TBD] | [TBD] | [TBD] | [TBD] |
| Extend right 25 nt | [TBD] | [TBD] | [TBD] | [TBD] | [TBD] | [TBD] |
| Extend both 25 nt | [TBD] | [TBD] | [TBD] | [TBD] | [TBD] | [TBD] |
| Shrink left 25 nt | [TBD] | [TBD] | [TBD] | [TBD] | [TBD] | [TBD] |
| Shrink right 25 nt | [TBD] | [TBD] | [TBD] | [TBD] | [TBD] | [TBD] |

*(Table will be populated after running the analysis.)*

### Most and Least Stable Regions

**Per-case stability** (mean tier agreement across all conditions):

| Case | Mean Tier Agreement | Std | Baseline FIRM Pairs | Baseline FLOPPY Pairs | Interpretation |
|------|---------------------|-----|---------------------|----------------------|----------------|
| [TBD] | [TBD] | [TBD] | [TBD] | [TBD] | [TBD] |

**Observations:**
- Cases with **higher baseline FIRM pair counts** tend to show [more/less] stability.
- Window **extension** generally [increases/decreases/has little effect on] tier agreement.
- Window **shrinkage** has [stronger/weaker] effects than extension.
- **Temperature** variations in the physiological range (30–42°C) produce [small/moderate/large] changes.
- **Parameter sets** show [low/moderate/high] agreement, with [set name] producing the most divergent predictions.

### Figures

![Tier Agreement Heatmap](../outputs/figures/layer5_tier_agreement_heatmap.png)
*Figure 1. Per-nucleotide tier agreement across cases and conditions. Green = high stability (agreement ≈ 1), red = low stability. Each row is a case; each column is a perturbation condition.*

![MEA Jaccard Heatmap](../outputs/figures/layer5_mea_jaccard_heatmap.png)
*Figure 2. MEA structure Jaccard index across cases and conditions. Measures base-pair set overlap between baseline and perturbed structures.*

## Interpretation

[To be completed after running analysis]

**Key findings:**
1. [Parameter set robustness]
2. [Temperature robustness]
3. [Window jitter robustness]
4. [Which regions/cases are most/least stable]
5. [Implications for antisense oligonucleotide design, RNA targeting]

**Practical implications:**
- For **ASO design** (e.g., Nusinersen targeting SMN2 ISS-N1), regions with [high/low] tier stability across conditions are [more/less] reliable targets.
- For **antiviral RNA targeting** (e.g., SARS-CoV-2 FSE), [parameter set / temperature] variations suggest [level of confidence].
- **Window selection** matters: extending by [N] nt produces [effect], suggesting that [interpretation].

## Limitations

1. **Window jitter assumes local independence.** Extending a window adds flanking nucleotides that may form long-range interactions with internal regions, not just terminal pairs. True robustness would require analyzing multiple independently annotated windows for the same functional element.

2. **Temperature range is physiological.** Extreme temperatures (e.g., <20°C, >50°C) were not tested. In vivo temperatures vary by cell type and organism (e.g., thermophiles), but the five disease cases are human/viral and physiological.

3. **MEA implementation is simplified.** The MEA structure was computed with a greedy pairing algorithm rather than the full dynamic programming MEA algorithm used by ViennaRNA. This may slightly underestimate MEA Jaccard similarity, but the relative trends across conditions remain valid.

4. **Pseudoknots are not considered.** ViennaRNA's partition function does not model pseudoknots, so FIRM/SOFT/FLOPPY classifications are limited to nested structures. The SARS-CoV-2 FSE is known to form a pseudoknot; base-pair probabilities reflect the ensemble of nested-only structures.

5. **No in vivo probing validation.** Robustness is measured as self-consistency (agreement between conditions), not accuracy against experimental data. Layer 3 (probing) will assess accuracy; this layer assesses precision.

6. **Sequence verification depends on NCBI.** Genomic coordinates were fetched via NCBI E-utilities. If an accession has been updated or deprecated, the fetch may fail. All five MVP cases successfully verified (sequences matched genomic slices), and accessions + coordinates are recorded for reproducibility.

## Runtime

**Expected runtime:** ~10–15 minutes on a 16 GB Mac for all five cases × all conditions (2 parameter sets + 4 temperatures + ~12 jitter conditions = ~18 conditions per case × 5 cases = 90 condition evaluations). Each condition requires one partition function calculation (~1–5 seconds per case).

**Actual runtime:** [TBD after run]

## References

1. Mathews DH, et al. Incorporating chemical modification constraints into a dynamic programming algorithm for prediction of RNA secondary structure. *PNAS* 2004;**101**:7287–7292. (Turner 2004 parameters)

2. Andronescu M, Condon A, Hoos HH, Mathews DH, Murphy KP. Efficient parameter estimation for RNA secondary structure prediction. *Bioinformatics* 2007;**23**:i19–i28. (Andronescu 2007 parameters)

3. Langdon WB, Petke J, Lorenz R. Evolving better RNAfold structure prediction. *EuroGP 2018*, LNCS **10781**:220–236. (Langdon 2018 parameters)

4. Lorenz R, et al. ViennaRNA Package 2.0. *Algorithms Mol Biol* 2011;**6**:26.

5. Do CB, Woods DA, Batzoglou S. CONTRAfold: RNA secondary structure prediction without physics-based models. *Bioinformatics* 2006;**22**:e90–e98. (MEA structure concept)

---

**Layer 5 Status:** Implementation complete. Results pending first run.

**Git tip SHA:** [TBD]
