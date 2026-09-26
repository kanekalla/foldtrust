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
| `sars2-fse` | SARS-CoV-2 frameshift element | COVID-19 / antiviral RNA | [DOI:10.1126/science.abf3546](https://doi.org/10.1126/science.abf3546) |
| `smn2-iss-n1` | SMN2 exon 7 and ISS-N1 | SMA / nusinersen target | [DOI:10.1056/NEJMoa1702752](https://doi.org/10.1056/NEJMoa1702752) |
| `cftr-5utr` | CFTR 5′ UTR and start codon region | Cystic fibrosis | [DOI:10.1126/science.2475911](https://doi.org/10.1126/science.2475911) |
| `mapt-e10` | MAPT exon 10 5′ss stem-loop | FTDP-17 / tauopathy | [DOI:10.1074/jbc.274.21.15134](https://doi.org/10.1074/jbc.274.21.15134) |
| `hcv-ires-dii` | HCV IRES domain II | Hepatitis C / viral IRES | [DOI:10.1128/JVI.73.2.1165-1174.1999](https://doi.org/10.1128/JVI.73.2.1165-1174.1999) |

Numbers below are from the committed demo outputs on `main` (`examples/out/*/report.md`). Re-running `foldtrust demo` after ViennaRNA upgrades can shift energies slightly; treat these as the recorded MVP snapshot.

---

## 5. Results

### 5.1 Summary table

**Note:** The original SMN2, CFTR, and MAPT sequences were invented (no NCBI provenance). HCV was the wrong domain (domain III, not II). All cases were rebuilt from exact NCBI coordinates (NG_008728.1:31999-32152, NM_000492.4:1-200, NG_007398.2:120818-121000, AF009606.1:44-118). See `docs/benchmark/layer0_cases.md` for provenance details.

| Case | Length | GC% | MFE (kcal/mol) | Firm | Soft | Floppy | Verdict |
|------|--------|-----|----------------|------|------|--------|---------|
| `smn2-iss-n1` | 154 | 30.5 | −31.30 | 2 | 6 | 1 | REDESIGN |
| `cftr-5utr` | 200 | 52.5 | −60.60 | 3 | 1 | 6 | REDESIGN |
| `mapt-e10` | 183 | 48.1 | −52.30 | 5 | 2 | 4 | REDESIGN |
| `hcv-ires-dii` | 75 | 54.7 | −23.80 | 2 | 3 | 0 | NEED PROBING |
| `sars2-fse` | 81 | 53.1 | −26.00 | 1 | 2 | 1 | REDESIGN |

Four of five cases land on **REDESIGN** (contain floppy stems). HCV IRES domain II has no floppy stems but soft stems outnumber firm, yielding **NEED PROBING**. These results demonstrate that ensemble reliability varies across disease windows: not all RNA structures are equally trustworthy under nearest-neighbor thermodynamics.

### 5.2 Case-by-case results

**Provenance note:** The original SMN2, CFTR, and MAPT sequences had no NCBI provenance (invented sequences). HCV was AF009606.1:148-415 (domain III into the CDS), not domain II (5'NTR 44-118). All were rebuilt from correct NCBI coordinates. Full audit and verification: `docs/benchmark/layer0_cases.md`.

#### SMN2 ISS-N1 (`smn2-iss-n1`) — NG_008728.1:31999-32152(+), 154 nt

- **MFE:** −31.30 kcal/mol; GC 30.5%; 9 stems.
- **Ensemble mix:** 2 FIRM, 6 SOFT, 1 FLOPPY.
- **Verdict:** REDESIGN (contains floppy stems).
- **Result in one line:** Local reliability is heterogeneous — exactly the situation ASO designers (nusinersen targets ISS-N1 at intron 7 +10..+27) should not flatten into one MFE cartoon.

#### CFTR 5′ UTR (`cftr-5utr`) — NM_000492.4:1-200(+), 200 nt

- **MFE:** −60.60 kcal/mol; GC 52.5%; 10 stems.
- **Ensemble mix:** 3 FIRM, 1 SOFT, 6 FLOPPY.
- **Verdict:** REDESIGN (contains floppy stems).
- **Result in one line:** A 5′ UTR with islands of firm structure but many floppy MFE helices — start codon accessibility is not reliably predicted from MFE alone.

#### MAPT exon 10 region (`mapt-e10`) — NG_007398.2:120818-121000(+), 183 nt

- **MFE:** −52.30 kcal/mol; GC 48.1%; 11 stems.
- **Ensemble mix:** 5 FIRM, 2 SOFT, 4 FLOPPY.
- **Verdict:** REDESIGN (contains floppy stems).
- **Result in one line:** The exon 10 / intron 10 5′ss stem-loop (involved in FTDP-17 pathogenic mutations) shows mixed reliability — firm stems coexist with floppy regions.

#### HCV IRES Domain II (`hcv-ires-dii`) — AF009606.1:44-118(+), 75 nt

- **MFE:** −23.80 kcal/mol; GC 54.7%; 5 stems.
- **Ensemble mix:** 2 FIRM, 3 SOFT, 0 FLOPPY.
- **Verdict:** NEED PROBING (no floppy stems, but soft stems outnumber firm).
- **Result in one line:** A structured viral IRES domain where most MFE stems have moderate-to-high ensemble support — useful contrast to the frameshift element.

#### SARS-CoV-2 frameshift element (`sars2-fse`) — NC_045512.2:13462-13542(+), 81 nt

- **MFE:** −26.00 kcal/mol; GC 53.1%; 4 stems.
- **Ensemble mix:** 1 FIRM, 2 SOFT, 1 FLOPPY.
- **Verdict:** REDESIGN (contains floppy stems).
- **Result in one line:** The programmed −1 ribosomal frameshift pseudoknot shows competing folds — trusting the MFE cartoon alone for antiviral design would be a mistake.

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
- SARS-CoV-2 frameshift: Kelly JA et al. *J Biol Chem* (2020). [DOI:10.1074/jbc.AC120.013449](https://doi.org/10.1074/jbc.AC120.013449); Bhatt PR et al. *Science* (2021). [DOI:10.1126/science.abf3546](https://doi.org/10.1126/science.abf3546)
- SMN2 ISS-N1: Singh NK et al. *Mol Cell Biol* (2006). [DOI:10.1128/MCB.26.4.1333-1346.2006](https://doi.org/10.1128/MCB.26.4.1333-1346.2006); nusinersen: Finkel RS et al. *NEJM* (2017). [DOI:10.1056/NEJMoa1702752](https://doi.org/10.1056/NEJMoa1702752)
- CFTR: Riordan JR et al. *Science* (1989). [DOI:10.1126/science.2475911](https://doi.org/10.1126/science.2475911); Zielenski J et al. *Genomics* (1991). [DOI:10.1016/0888-7543(91)90503-7](https://doi.org/10.1016/0888-7543(91)90503-7)
- MAPT exon 10: Grover A et al. *J Biol Chem* (1999). [DOI:10.1074/jbc.274.21.15134](https://doi.org/10.1074/jbc.274.21.15134); Varani L et al. *PNAS* (1999). [DOI:10.1073/pnas.96.14.8229](https://doi.org/10.1073/pnas.96.14.8229); Hutton M et al. *Nature* (1998). [DOI:10.1038/31508](https://doi.org/10.1038/31508)
- HCV IRES domain II: Honda M et al. *J Virol* (1999). [DOI:10.1128/JVI.73.2.1165-1174.1999](https://doi.org/10.1128/JVI.73.2.1165-1174.1999); Lukavsky PJ et al. *Nat Struct Biol* (2003). [DOI:10.1038/nsb1004](https://doi.org/10.1038/nsb1004)

Full case-level provenance (NCBI accessions, coordinates, landmark tables) in `data/cases/*/meta.yaml` and `docs/benchmark/layer0_cases.md`.

---

## 9. Benchmark Results Summary

**Updated:** September 25, 2026

FoldTrust's six benchmark layers summarise what the six benchmark layers show about per-stem confidence tiers (FIRM/SOFT/FLOPPY). Full results and methods in `BENCHMARK.md`.

### Key Findings

**Layer 1 (Scoring correctness):** 5/5 scoring tests pass (perfect/no/half overlap, ECE, FIRM tier accuracy).

**Layer 2 (Structure accuracy):** MEA F1 = 0.563 (95% CI: 0.542–0.585) on 600 ArchiveII/Rfam/bpRNA structures. Moderate agreement with comparative structures reflects thermodynamic vs. phylogenetic modeling differences.

**Layer 3 (Calibration):** ECE = 0.0664; AUROC = 0.8889. Tier PPV: FIRM 0.674, SOFT 0.299, FLOPPY 0.146 (pooled over MFE pairs). Higher tiers had higher reference agreement (PPV FIRM 0.674 > SOFT 0.299 > FLOPPY 0.146).

**Layer 4 (SHAPE agreement):** Spearman ρ = 0.29–0.55 for 4 of 5 datasets (Pyle 0.16, n.s.); genome control 43rd–82nd percentile, empirical p 0.18–0.57 (not exceptional).

**Layer 5 (Robustness):** FIRM tier retention: 95% at 25/30°C, 100% at 42°C (pooled over 5 disease cases). Andronescu2007 80%, Langdon2018 60% (vs. Turner2004 baseline). Window context: 50-100 nt flanks partially recover structure (MEA BP distance 19.6–22.4).

### Interpretation

Higher tiers had higher reference agreement within ViennaRNA's ensemble predictions. FIRM stems (P ≥ 0.85) have 67% PPV against reference structures and robust retention across temperatures. Lower layers (SOFT/FLOPPY) are increasingly uncertain. The system is a **hypothesis generator** for which predicted stems to prioritize; experimental validation remains essential.

**Limitations:** Thermodynamic model only (no phylogeny, no RBPs); SHAPE coverage limited to FSE; FLOPPY tier uninformative (PPV 0.15); parameter sensitivity (Andronescu/Langdon shift tiers); context dependence (flanks alter predictions).

### Reproduce

```bash
# Fetch data bundle (see MANIFEST.md for URLs or extract from uploads/)
mkdir -p data/_cache
tar -xzf uploads/foldtrust_bench_data.tar.gz -C data/_cache
tar -xzf uploads/bprna_TS0_canonicals.tar.gz -C data/_cache

# Run all layers (~3 minutes, 16GB machine)
foldtrust benchmark all --output benchmarks/outputs

# Individual layers
foldtrust benchmark layer1  # Unit tests
foldtrust benchmark layer2  # Structure accuracy
foldtrust benchmark layer3  # Calibration
foldtrust benchmark layer4_shape  # SHAPE agreement
foldtrust benchmark layer5  # Robustness
```

All outputs saved to CSV/JSON in `benchmarks/outputs/layerN/`. Figures in `benchmarks/outputs/figures/`.

**Drug-discovery relevance:** See `docs/benchmark/impact.md` for hypothesis-level discussion of ASO/small-molecule targeting (SMN2 ISS-N1, SARS-CoV-2 FSE, HCV IRES, CFTR 5'UTR, MAPT exon 10). No clinical claims; validation required.

**Data sources:** ArchiveII (Sloma & Mathews 2016), Rfam 15.1, bpRNA-1m, SARS-CoV-2 SHAPE (DasLab repo). See `MANIFEST.md` and `SHA256SUMS`.

**Code:** `src/foldtrust/benchmark/`, tests in `tests/`, docs in `docs/benchmark/`.

---

## 10. References (methods + cases + benchmark)
