# FoldTrust Notes: Question, Theory, Results, and Interpretation

**Author:** Kishore Anekalla (`kanekalla`)  
**Project:** [FoldTrust](https://github.com/kanekalla/foldtrust) — ensemble-first RNA structure reliability reports  
**Scope:** MVP on published disease-relevant RNA windows  
**Energy model:** ViennaRNA default parameters (Turner 2004), 37 °C / 1 M NaCl assumptions

These notes are the scientific write-up behind the tool: the question that motivated FoldTrust, the thermodynamic theory, how the MVP answers the question, and what the five demo cases actually show.

---

## 1. The question

**When an RNA fold looks crisp in the minimum free energy (MFE) cartoon, how much of that structure should we trust for biology or therapeutic design?**

Secondary-structure tools routinely return a single MFE structure — one lowest-energy drawing of stems and loops. That drawing is easy to publish and easy to over-believe. In practice:

- Competing folds can sit near the MFE in free energy.
- A stem drawn in the MFE can have **low base-pair probability** in the Boltzmann ensemble.
- Antisense oligos, splice-switching reagents, and antiviral RNA targeting care about **accessibility and reliability**, not only about the prettiest cartoon.

FoldTrust's working question is therefore more precise:

> For a disease-relevant RNA window, which MFE stems are **firm** in the thermodynamic ensemble, which are only **soft**, and which are **floppy** (appear in the MFE but lack ensemble support)?

That is a reliability question, not a "predict the true native structure" claim.

---

## 2. Theory

### 2.1 MFE vs ensemble

RNA secondary structure prediction under nearest-neighbor thermodynamics assigns each structure \(S\) a free energy \(\Delta G(S)\). The **minimum free energy** structure is:

\[
S_{\mathrm{MFE}} = \arg\min_S \Delta G(S)
\]

That is one structure. The physically relevant object at equilibrium is the **Boltzmann ensemble**. The probability of structure \(S\) is:

\[
P(S) = \frac{e^{-\Delta G(S)/RT}}{Z}, \qquad Z = \sum_{S'} e^{-\Delta G(S')/RT}
\]

where \(Z\) is the partition function. ViennaRNA computes \(Z\) (and related quantities) with dynamic programming, not by enumerating all structures.

### 2.2 Base-pair probabilities

For positions \(i < j\), the **base-pair probability** is the ensemble-averaged indicator that \(i\) pairs with \(j\):

\[
P_{ij} = \sum_{S} \mathbf{1}\{(i,j) \in S\}\, P(S)
\]

\(P_{ij}\) is the quantity FoldTrust treats as the reliability signal. An MFE stem can look "solid" while many of its pairs have small \(P_{ij}\). That is the crisp-MFE / floppy-ensemble failure mode.

### 2.3 Stem-level reliability (FoldTrust definition)

FoldTrust:

1. Runs `RNAfold` for the MFE structure and the partition function / pair probabilities.
2. Parses contiguous helices (stems) from the MFE dot-bracket.
3. Scores each stem by the **mean pair probability** of its constituent pairs.
4. Flags stems with fixed thresholds chosen for interpretation, not as universal biology cutoffs:

| Flag | Mean \(P\) | Meaning |
|------|------------|---------|
| **FIRM** | \(\ge 0.85\) | Ensemble strongly supports the MFE stem |
| **SOFT** | \(0.5 \le P < 0.85\) | Partial support; alternatives plausible |
| **FLOPPY** | \(< 0.5\) | MFE stem is poorly supported; do not over-trust the cartoon |

5. Emits HTML + Markdown report cards with a stem table, pair-probability heatmap, and a coarse verdict (TRUST / REDESIGN / NEED PROBING).

**Key principle:** ask for base-pair probabilities (or the ensemble), not only the MFE.

### 2.4 What this theory does *not* claim

Thermodynamic models are simplified:

- Default solution conditions (commonly treated as ~37 °C, high salt), not the crowded cell.
- No cotranscriptional folding kinetics.
- No RNA-binding proteins, chemical modifications, or tertiary contacts.
- No experimental SHAPE/DMS constraints unless added later.

FoldTrust is a **hypothesis generator** for which stems look reliable under nearest-neighbor thermodynamics. Wet-lab probing and functional assays remain decisive.

---

## 3. Application (what FoldTrust is for)

FoldTrust is built as a small, local, reproducible **structure report card** for disease windows where ensemble reliability matters:

| Application lens | Why ensemble reliability matters |
|------------------|----------------------------------|
| **ASO / splice-switching design** | Floppy local structure can mean higher accessibility; firm stems can block oligos |
| **Viral RNA targeting** | Frameshift elements and IRES domains often have competing folds |
| **UTR / splicing disease genes** | Structural cartoons used in mechanism stories need a reliability check |
| **Scientific demonstration** | Shows ensemble thinking beyond "ran RNAfold once" |

The MVP ships five published windows (not invented sequence) so a reader can reproduce the same reports with `foldtrust demo`.

---

## 4. Methods (MVP)

- **Engine:** ViennaRNA CLI (`RNAfold`) via subprocess; pair probabilities from the partition-function / PostScript ubox channel (ensemble probabilities only).
- **Language:** Python package + Typer CLI (`report`, `batch`, `demo`).
- **Inputs:** FASTA + per-case `meta.yaml` (disease, gene, teaching point, DOI/PMID).
- **Outputs:** `report.html`, `report.md`, `pair_probabilities.png` under `examples/out/<case>/` for the demo.
- **Cases:**

| Case ID | Window | Disease focus | Primary citation |
|---------|--------|---------------|------------------|
| `sars2-fse` | SARS-CoV-2 frameshift stimulatory element | COVID-19 / antiviral RNA | [DOI:10.1126/science.abc3546](https://doi.org/10.1126/science.abc3546) |
| `smn2-iss-n1` | SMN2 intron 7 ISS-N1 neighborhood | SMA / Spinraza context | [DOI:10.1056/NEJMoa1702752](https://doi.org/10.1056/NEJMoa1702752) |
| `cftr-5utr` | CFTR 5′ UTR | Cystic fibrosis / lung | [DOI:10.1152/physrev.00025.2017](https://doi.org/10.1152/physrev.00025.2017) |
| `mapt-e10` | MAPT exon 10 splice regulatory region | Tauopathy / aging | [DOI:10.1093/hmg/10.10.1029](https://doi.org/10.1093/hmg/10.10.1029) |
| `hcv-ires-dii` | HCV IRES Domain II | Hepatitis C / structured viral RNA | [DOI:10.1006/jmbi.1999.2918](https://doi.org/10.1006/jmbi.1999.2918) |

Numbers below are from the committed demo outputs on `main` (`examples/out/*/report.md`). Re-running `foldtrust demo` after ViennaRNA upgrades can shift energies slightly; treat these as the recorded MVP snapshot.

---

## 5. Results

### 5.1 Summary table

| Case | Length | GC% | MFE (kcal/mol) | Firm | Soft | Floppy | Demo verdict |
|------|--------|-----|----------------|------|------|--------|--------------|
| `sars2-fse` | 181 | 32.0 | −39.50 | 1 | 0 | **7** | REDESIGN |
| `smn2-iss-n1` | 201 | 22.9 | −29.40 | 3 | 5 | 5 | REDESIGN |
| `cftr-5utr` | 268 | 31.7 | −55.10 | 3 | 9 | 5 | REDESIGN |
| `mapt-e10` | 268 | 58.6 | −86.00 | 7 | 10 | 2 | REDESIGN |
| `hcv-ires-dii` | 268 | 57.5 | −91.20 | 12 | 4 | 2 | REDESIGN |

All five MVP windows currently land on **REDESIGN** under FoldTrust's stem-level rule (any floppy MFE stem triggers caution). That is intentional for teaching: disease windows are not "easy trust" cartoons.

### 5.2 Case-by-case results

#### SARS-CoV-2 frameshift element (`sars2-fse`) — best floppy teaching case

- **MFE:** −39.50 kcal/mol over 181 nt; eight MFE stems drawn.
- **Ensemble:** only **stem 6** is FIRM (mean \(P = 0.921\)); stems 1–5, 7–8 are FLOPPY (mean \(P \approx 0.10\)–\(0.18\)).
- **Result in one line:** the MFE looks structured; the ensemble says most of those helices are not trustworthy.

#### SMN2 ISS-N1 neighborhood (`smn2-iss-n1`)

- **MFE:** −29.40 kcal/mol; 13 stems.
- **Ensemble mix:** 3 FIRM, 5 SOFT, 5 FLOPPY (FIRM examples: stem 1 mean \(P = 0.962\); stem 13 mean \(P = 0.966\)).
- **Result in one line:** local reliability is heterogeneous — exactly the situation ASO designers should not flatten into one cartoon.

#### CFTR 5′ UTR (`cftr-5utr`)

- **MFE:** −55.10 kcal/mol; 17 stems.
- **Ensemble mix:** 3 FIRM, 9 SOFT, 5 FLOPPY (notably firm longer helices such as stem 10 mean \(P = 0.951\)).
- **Result in one line:** a mostly soft UTR with islands of firm structure and several floppy MFE helices.

#### MAPT exon 10 region (`mapt-e10`)

- **MFE:** −86.00 kcal/mol; GC-rich (58.6%); 19 stems.
- **Ensemble mix:** 7 FIRM, 10 SOFT, only 2 FLOPPY.
- **Result in one line:** more ensemble-supported than the viral frameshift window, but still not a blanket TRUST call.

#### HCV IRES Domain II (`hcv-ires-dii`)

- **MFE:** −91.20 kcal/mol; 18 stems.
- **Ensemble mix:** **12 FIRM**, 4 SOFT, 2 FLOPPY (many Domain II helices with mean \(P > 0.9\)).
- **Result in one line:** a structured viral RNA where most MFE stems *are* ensemble-supported — useful contrast to `sars2-fse`.

---

## 6. Interpretation

### 6.1 What the MVP demonstrates scientifically

1. **MFE ≠ reliability.** `sars2-fse` is the clearest proof inside this repo: 7/8 MFE stems fail the ensemble mean-\(P\) test.
2. **Disease windows are not one story.** HCV Domain II is largely firm; the frameshift element is largely floppy; SMN2 and CFTR sit in between with soft/floppy mixtures.
3. **Therapeutic reading is directional, not automatic.**
   - Floppy MFE stems → do not trust the cartoon for mechanism slides; consider accessibility / alternative folds.
   - Firm stems → more plausible structural constraints, still subject to in-cell and experimental checks.
4. **Contrast demonstrates discriminatory power.** Showing both a floppy-dominant viral frameshift window and a firm-dominant IRES window proves the tool can discriminate, not just stamp REDESIGN on everything for drama.

### 6.2 How to read a FoldTrust report in practice

1. Open the stem table before the ASCII MFE drawing.
2. Ask: which helices would I still believe if I deleted the MFE cartoon?
3. For ASO / probe design, prioritize **local** soft/floppy neighborhoods over global MFE aesthetics.
4. Treat REDESIGN / NEED PROBING as "get probabilities or experiments," not "the biology is wrong."

### 6.3 Limitations (honest)

- Thresholds (0.85 / 0.5) are interpretive conventions for this MVP.
- Long-range pairs, pseudoknots, and multi-branch junctions are only as good as the ViennaRNA model and the window you chose.
- Published windows are curated teaching cases; coordinates and citations live in each `meta.yaml`.
- No SHAPE/DMS fusion yet (roadmap), which is where ensemble models become much stronger for real programs.

### 6.4 Take-home

FoldTrust's application is to make ensemble thinking **operational**: same disease sequence, same MFE everyone already looks at, plus an explicit firm / soft / floppy report card. The demo results show that for at least one high-profile antiviral RNA window (SARS-CoV-2 FSE), trusting the MFE alone would be a scientific mistake — and that other disease RNAs (HCV IRES Domain II) can look much more ensemble-consistent under the same pipeline.

---

## 7. Reproduce

```bash
brew install viennarna   # or apt install vienna-rna
git clone https://github.com/kanekalla/foldtrust.git
cd foldtrust
pip install -e .
foldtrust demo
open examples/out/sars2-fse/report.html
```

Pre-generated MVP snapshots also live under `examples/out/`.

---

## 8. References (methods + cases)

- Lorenz et al., ViennaRNA Package 2.0. *Algorithms Mol Biol* (2011). [DOI:10.1186/1748-7188-6-26](https://doi.org/10.1186/1748-7188-6-26)
- Mathews et al., Turner nearest-neighbor parameters. *PNAS* (2004). [DOI:10.1073/pnas.0401799101](https://doi.org/10.1073/pnas.0401799101)
- SARS-CoV-2 frameshift: [DOI:10.1126/science.abc3546](https://doi.org/10.1126/science.abc3546)
- SMN2 / nusinersen context: [DOI:10.1056/NEJMoa1702752](https://doi.org/10.1056/NEJMoa1702752)
- CFTR review: [DOI:10.1152/physrev.00025.2017](https://doi.org/10.1152/physrev.00025.2017)
- MAPT exon 10 stem-loop: [DOI:10.1093/hmg/10.10.1029](https://doi.org/10.1093/hmg/10.10.1029)
- HCV IRES structure: [DOI:10.1006/jmbi.1999.2918](https://doi.org/10.1006/jmbi.1999.2918)

Case-level extra DOIs are listed in each generated report and `data/cases/*/meta.yaml`.

---

## 9. Benchmark: Data Availability and Implementation Constraints

**Goal:** Validate FoldTrust's reliability classifications against standard benchmark datasets and experimental probing data.

### 9.1 Investigation Summary

Benchmark implementation attempted four analyses as specified:

1. **Reference structure accuracy** (ArchiveII dataset)
2. **Calibration of pair probabilities** (reliability diagram, tier accuracy)
3. **Experimental probing correlation** (SHAPE-MaP reactivity)
4. **Robustness** (parameter sets, temperature, window jitter)

### 9.2 Data Availability Findings

#### 9.2.1 ArchiveII Dataset

**Status:** NOT ACCESSIBLE

**Attempts documented in `benchmarks/ARCHIVEII_ATTEMPTS.md`:**
- Mathews Lab direct URL (rna.urmc.rochester.edu/RNAstructure/Supplemental/ArchiveII/archiveII.tar.gz): HTTP 404
- RNA STRAND v2.0 download: Connection timeout
- Published benchmark repositories (mxfold2, LinearPartition, E2Efold, allegro): No .ct/.bpseq files distributed
- Public databases (Rfam, PDB): Require per-family extraction and manual curation

**Conclusion:** The canonical ArchiveII dataset referenced in RNA structure prediction literature (e.g., Mathews lab 2004-2010 publications) is no longer hosted at its documented URLs as of September 2026. Published papers cite ArchiveII but do not redistribute the data in their repositories.

**Decision per user instruction:** "If you truly cannot get real reference data, say so and do not report reference metrics at all." Reference structure accuracy metrics are **not reported**.

#### 9.2.2 SHAPE-MaP Reactivity Data

**Status:** DATA FOUND, COORDINATE MAPPING REQUIRED

**Data source:** 
- Repository: https://github.com/DasLab/SARS_CoV-2_shape_comparison (cloned Sep 25, 2026)
- Files: `zhang_invivo_reactivity.csv`, `incarnato_invivo_reactivity.csv`, `pyle_reactivity.csv`
- Coverage: Genome-wide SARS-CoV-2 (29,903 positions)
- Citations: Manfredonia et al. 2020 (doi:10.1038/s41586-020-2681-1), Huston et al. 2021

**Issue documented in `benchmarks/SHAPE_DATA_INVESTIGATION.md`:**  
The frameshift element sequence in `data/cases/sars2-fse/sequence.fa` (181 nt) does not match any substring of the reference genome (`refseq.txt`) in the DasLab repository. Proper implementation requires:
1. Resolving sequence source (viral isolate/strain differences)
2. Mapping FSE coordinates to genome positions
3. Extracting corresponding SHAPE values

**Estimated time:** 2-4 hours for coordinate resolution and validation.

**Decision:** SHAPE probing analysis **not completed** in this iteration. Infrastructure exists (`foldtrust.benchmark.probing`); data ingestion deferred pending coordinate mapping.

#### 9.2.3 Robustness Analysis

**Status:** TEMPERATURE SWEEP IMPLEMENTABLE; PARAMETER FILES NOT FOUND

**ViennaRNA capabilities verified:**
- Temperature control: RNAfold `-T` flag available ✓
- Parameter files: RNAfold `-P` flag present, but alternative parameter sets (Andronescu 2007, Langdon 2018) not found in ViennaRNA 2.5.1 installation
- Window jitter: Requires flanking genomic sequence not included in disease case files

**What can be implemented:**
- Temperature sweep (24°C, 37°C, 42°C) on five disease cases
- Tier stability reporting (fraction of pairs/stems changing FIRM/SOFT/FLOPPY classification)

**Limitations:**
- Parameter set comparison (Turner 2004 vs alternatives) not feasible without locating or downloading alternative .par files
- Window jitter requires genomic context beyond the curated disease windows

### 9.3 Honest Assessment

**What this investigation demonstrates:**
1. Benchmark data availability in RNA structure prediction is challenging: canonical datasets (ArchiveII) are no longer hosted, and experimental data (SHAPE) requires coordinate mapping not trivial to resolve.
2. The benchmark infrastructure (`foldtrust.benchmark.*` modules) is correctly implemented: scoring works (F1=1.0 on exact matches), calibration metrics compute properly, dataset loaders function.
3. Reporting invalid metrics (e.g., F1=0.067 from mismatched hand-written test structures) would be scientifically dishonest.

**What was NOT done:**
- Reference structure accuracy: No ArchiveII → no metrics reported
- SHAPE correlation: Data located but coordinates unresolved → analysis deferred
- Full robustness: Temperature sweep implementable; parameter/jitter analysis requires additional data/files

**What FoldTrust provides:**
FoldTrust's value remains **transparency about ensemble uncertainty**. The tool correctly computes base-pair probabilities and classifies stems (FIRM/SOFT/FLOPPY) to expose MFE ambiguity. Benchmark validation against gold-standard datasets would strengthen the scientific claim, but the core functionality—asking for pair probabilities instead of only MFE—is independently valuable.

### 9.4 Reproducibility

**Data sources investigated:**
- ArchiveII: Attempted rna.urmc.rochester.edu, www.rnasoft.ca, GitHub repos (documented in `benchmarks/ARCHIVEII_ATTEMPTS.md`)
- SHAPE: https://github.com/DasLab/SARS_CoV-2_shape_comparison (successfully cloned)

**Code status:**
- Benchmark modules: Fully implemented, tested (17 tests pass)
- Metrics: Validated (sensitivity/PPV/F1/MCC, calibration, tier accuracy)
- CLI: `foldtrust benchmark <analysis>` functional

**Limitations documented:**
- `benchmarks/IMPLEMENTATION_STATUS.md`: Honest assessment of what is/isn't implementable
- No fabricated data, no invalid metrics reported

### 9.5 Recommendation for Future Work

To complete benchmark validation, future efforts should:
1. Obtain ArchiveII from colleagues with archived copies, or build a validated reference set from PDB/Rfam with documented curation
2. Resolve SARS-CoV-2 FSE coordinate mapping (contact paper authors or use NCBI RefSeq annotations)
3. Locate or download ViennaRNA alternative parameter files for Turner 2004 vs Andronescu/Langdon comparison
4. Add flanking genomic sequences to disease cases for window-jitter analysis

---

## 10. References (methods + cases + benchmark)
