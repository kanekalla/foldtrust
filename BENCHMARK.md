# FoldTrust Benchmark Report

**Date:** September 25, 2026  
**Runtime:** ~3 minutes for default subsets (16GB machine)

This document reports FoldTrust's six-layer benchmark for RNA secondary structure reliability classification. All numbers are computed from saved outputs in `benchmarks/outputs/` — no hand-typed results.

---

## What FoldTrust Claims

FoldTrust classifies predicted RNA base pairs into three reliability tiers based on ensemble probabilities from ViennaRNA:

- **FIRM** (P ≥ 0.85): High-confidence stems
- **SOFT** (0.5 ≤ P < 0.85): Moderate-confidence stems  
- **FLOPPY** (P < 0.5): Low-confidence stems

The benchmark tests these tiers using:
1. Unit tests on scoring logic correctness
2. Agreement with comparative/crystallographic reference structures
3. Calibration of predicted probabilities
4. Correlation with SHAPE chemical probing
5. Robustness to temperature, parameters, and window context

---

## Synthesis: Benchmark Layer Summary

| Layer | Question | Data (n) | Headline Metric(s) | Verdict | Doc Link |
|-------|----------|----------|-------------------|---------|----------|
| **0** | Are disease windows correctly defined? | 5 cases | FIRM/SOFT/FLOPPY stems per case | All coordinates verified from NCBI records | [layer0_cases.md](docs/benchmark/layer0_cases.md) |
| **1** | Is tier-calling logic correct? | 5 unit tests | 5/5 passed | Scoring logic passed 5/5 unit tests | [layer1_scoring.md](docs/benchmark/layer1_scoring.md) |
| **2** | Accuracy vs. reference structures? | 600 structures (200 default) | F1: MFE 0.548, MEA 0.563, Centroid 0.570 | Low absolute agreement reflects thermodynamic vs. comparative differences | [layer2_accuracy.md](docs/benchmark/layer2_accuracy.md) |
| **3** | Are probabilities calibrated? | 185,653 pairs (21,191 positive) | ECE 0.0664; FIRM PPV 0.674, SOFT 0.299, FLOPPY 0.146 | Higher tiers have higher reference agreement; high-prob pairs overconfident (ECE p≥0.5 = 0.317) | [layer3_calibration.md](docs/benchmark/layer3_calibration.md) |
| **4** | SHAPE correlation? | FSE + 5 SHAPE datasets | Spearman ρ 0.29–0.55 (4 of 5 significant); AUROC 0.64–0.86 | Weak-to-moderate correlation; Pyle ρ=0.16 n.s.; genome control 43rd–82nd percentile (p 0.18–0.57) | [layer4_shape.md](docs/benchmark/layer4_shape.md) |
| **5** | Robustness? | 5 cases × (3 temperatures + 3 parameter sets + 4 flank sizes) | FIRM retention: 0.95 at 25/30°C, 0.60–0.80 alternate params, 0.45–0.70 with flanks | FIRM stems mostly robust; context matters (25-nt flanks disrupt many) | [layer5_robustness.md](docs/benchmark/layer5_robustness.md) |

---

## Benchmark Layers

### Layer 0: Disease Case Definitions

**Question:** Are the five disease-relevant RNA windows correctly defined?

**Data:** 5 cases in `data/cases/` (see `docs/benchmark/layer0_cases.md`)  
**Source:** `benchmarks/outputs/layer0_cases/case_metrics.csv`

| Case | Coordinates | Length | FIRM / SOFT / FLOPPY | Verdict |
|------|-------------|--------|---------------------|---------|
| **sars2-fse** | NC_045512.2:13462-13542 (+) | 81 nt | 1 / 2 / 1 | REDESIGN |
| **smn2-iss-n1** | NG_008728.1:31999-32152 (+) | 154 nt | 2 / 6 / 1 | REDESIGN |
| **cftr-5utr** | NM_000492.4:1-200 (+) | 200 nt | 3 / 1 / 6 | REDESIGN |
| **mapt-e10** | NG_007398.2:120818-121000 (+) | 183 nt | 5 / 2 / 4 | REDESIGN |
| **hcv-ires-dii** | AF009606.1:44-118 (+) | 75 nt | 2 / 3 / 0 | NEED PROBING |

**Verdict:** Coordinates re-derived from NCBI records (RefSeq NC_/NG_/NM_; GenBank AF009606.1).

([Full report](docs/benchmark/layer0_cases.md))

---

### Layer 1: Scoring Correctness

**Question:** Does the tier-calling logic work correctly?

**Method:** Unit tests on synthetic sequences  
**Command:** `foldtrust benchmark layer1`  
**Source:** `benchmarks/outputs/layer1/layer1_tests.json`

**Results:** 5/5 scoring tests passed

| Test | Result |
|------|--------|
| Perfect match | Sensitivity 1.0, PPV 1.0, F1 1.0 ✓ |
| No overlap | Sensitivity 0.0, PPV 0.0, F1 0.0 ✓ |
| Half overlap | Sensitivity 0.5, PPV 0.5, F1 0.5 ✓ |
| Calibration (ECE) | ECE 0.0 ✓ |
| Tier accuracy (FIRM) | FIRM PPV 1.0 ✓ |

**Parameter-loading energies** are tested in `tests/test_param_energies.py` (FSE MFE = -26.0 Turner2004, -22.26 Andronescu2007, -24.70 Langdon2018).

**Verdict:** Tier logic is correct.

([Full tests](benchmarks/outputs/layer1/layer1_tests.json))

---

### Layer 2: Structure Accuracy

**Question:** How well do FoldTrust's MFE/MEA/centroid structures agree with comparative/crystallographic references?

**Data:** 600 structures from ArchiveII, Rfam, bpRNA (200-structure default subset)  
**Source:** `data/_cache/` (verified SHA256, see `MANIFEST.md`)  
**Command:** `foldtrust benchmark layer2`

**Results (exact base-pair matching; mean of per-structure metrics; `benchmarks/outputs/layer2/layer2_summary_overall.csv`):**

| Method | Sensitivity | PPV | F1 | MCC |
|--------|-------------|-----|----|----|
| **MFE** | 0.639 (0.616–0.662) | 0.499 (0.477–0.520) | 0.548 (0.526–0.569) | 0.557 (0.536–0.578) |
| **MEA** | 0.641 (0.619–0.662) | 0.521 (0.500–0.543) | **0.563** (0.542–0.585) | 0.571 (0.550–0.593) |
| **Centroid** | 0.622 (0.600–0.643) | 0.551 (0.530–0.574) | 0.570 (0.549–0.590) | 0.577 (0.556–0.598) |

(95% bootstrap CI, 1000 iterations)

**1-nt slip tolerance:** F1 0.579/0.594/0.599 (MFE/MEA/Centroid), MCC 0.589/0.603/0.607 (`layer2_raw_results.csv`)

**Verdict:** MEA/centroid slightly better F1 than MFE (~0.56–0.57). Low absolute agreement reflects known thermodynamic vs. comparative structure differences.

([Full results](benchmarks/outputs/layer2/), [Detailed docs](docs/benchmark/layer2_accuracy.md))

---

### Layer 3: Calibration

**Question:** Are the predicted pair probabilities well-calibrated? Do FIRM/SOFT/FLOPPY tiers correctly rank reliability?

**Data:** Same 600 structures as Layer 2  
**Method:** Expected Calibration Error (ECE), AUROC, tier-level PPV  
**Command:** `foldtrust benchmark layer3`

**Results:**

| Metric | Value | Interpretation |
|--------|-------|----------------|
| **ECE** | 0.0664 (p≥0.5: 0.3174) | Probabilities moderately calibrated; high-prob pairs overconfident |
| **AUROC** | 0.8889 | Good discrimination of true vs. false pairs |
| **AUPRC** | 0.6153 | Precision-recall area; reflects imbalance |

**Tier PPV (MFE pairs):**

| Tier | Pooled PPV | Structures | Total Pairs | Interpretation |
|------|-----------|------------|-------------|----------------|
| **FIRM** | **0.674** | 582 | 16,417 | FIRM pairs match the reference 67% of the time |
| **SOFT** | **0.299** | 585 | 8,518 | SOFT pairs match 30% of the time |
| **FLOPPY** | **0.146** | 466 | 5,504 | FLOPPY pairs match 15% of the time |

**Verdict:** Higher tiers had higher reference agreement (PPV FIRM 0.674 > SOFT 0.299 > FLOPPY 0.146). High-probability pairs are over-confident (ECE p≥0.5 = 0.317). See calibration plots in `benchmarks/outputs/figures/`.

([Full results](benchmarks/outputs/layer3/layer3_summary.json), [Docs](docs/benchmark/layer3_calibration.md))

---

### Layer 4: SHAPE Agreement

**Question:** Do unpaired probabilities correlate with SHAPE reactivity on SARS-CoV-2 FSE?

**Data:** 5 SHAPE/icSHAPE datasets (Incarnato, Pyle, Zhang) on NC_045512.2 (29,903 nt)  
**Method:** Spearman correlation between unpaired probability and normalized reactivity; AUROC for reactive vs. unreactive classification  
**Sequence:** FSE NC_045512.2:13462-13542 (81 nt)  
**Command:** `foldtrust benchmark layer4_shape`

**Results:**

| Dataset | Chemistry | Spearman ρ | p-value | AUROC | Coverage |
|---------|-----------|------------|---------|-------|----------|
| **incarnato_invitro** | SHAPE-MaP | **0.294** | 0.0077 | 0.705 | 100% |
| **incarnato_invivo** | SHAPE-MaP | **0.355** | 0.0013 | 0.785 | 98.8% |
| **pyle** | SHAPE-MaP | 0.163 | 0.145 | 0.644 | 100% |
| **zhang_invitro** | icSHAPE | **0.497** | 2.4e-06 | 0.850 | 100% |
| **zhang_invivo** | icSHAPE | **0.547** | 1.2e-07 | 0.864 | 100% |

**Genome control:** FSE correlation is not exceptional: 43rd–82nd percentile of 363–369 non-overlapping 81-nt windows, empirical p 0.18–0.57 (`layer4_shape_results.json`).

**Verdict:** Weak-to-moderate correlation for 4 of 5 datasets (ρ 0.29–0.55; Pyle ρ 0.16, n.s.); AUROC 0.64–0.86. SHAPE-directed folding: MEA F1 0.83–0.98, MFE F1 0.40–0.96 vs unconstrained (`shape_directed_folding.csv`).

([Full results](benchmarks/outputs/layer4_shape/), [Figures](benchmarks/outputs/figures/layer4*.png), [Docs](docs/benchmark/layer4_shape.md))

---

### Layer 5: Robustness

**Question:** Are tier classifications robust to temperature, parameter sets, and window context?

**Data:** 5 disease cases (FSE, SMN2, CFTR, MAPT, HCV)  
**Method:** Compare tier assignments at 25/30/42°C vs. 37°C baseline; Andronescu2007/Langdon2018 vs. Turner2004; flanking context 0/25/50/100 nt  
**Command:** `foldtrust benchmark layer5`

**Results:**

#### Temperature Sweep (pooled retention = retained / total reference stems)

| Tier | 37°C Baseline | 25°C | 30°C | 42°C |
|------|--------------|------|------|------|
| **FIRM** | 1.000 | 0.950 | 0.950 | 1.000 |
| **SOFT** | 0.733 | 0.600 | 0.600 | 0.733 |
| **FLOPPY** | 0.364 | 0.273 | 0.364 | 0.364 |

#### Parameter Set Sweep

| Tier | Turner2004 Baseline | Andronescu2007 | Langdon2018 |
|------|---------------------|----------------|-------------|
| **FIRM** | 1.000 | 0.800 | 0.600 |
| **SOFT** | 0.733 | 0.067 | 0.133 |
| **FLOPPY** | 0.364 | 0.455 | 0.091 |

#### Window Context (flanking nucleotides)

| Tier | 0 nt | 25 nt | 50 nt | 100 nt |
|------|------|-------|-------|--------|
| **FIRM** | 1.000 | 0.450 | 0.700 | 0.700 |
| **SOFT** | 0.733 | 0.200 | 0.200 | 0.133 |
| **FLOPPY** | 0.364 | 0.545 | 0.364 | 0.455 |

**Verdict:** FIRM tier is robust to temperature (95% retention at 25/30°C) but sensitive to alternative parameter sets (60-80% retention). Context matters: 25-nt flanks disrupt many stems (MEA BP distance = 30.4); 50-100 nt partially recover (distance = 19.6–22.4). SOFT/FLOPPY tiers are less stable.

([Full results](benchmarks/outputs/layer5/), [Figures](benchmarks/outputs/figures/layer5*.png), [Docs](docs/benchmark/layer5_robustness.md))

---

## Disease Window Summary

| Case | Tier Summary (case_metrics.csv) | Layer 4 | Layer 5 | Caveats | Link |
|------|----------------------------------|---------|---------|---------|------|
| **SARS-CoV-2 FSE** | 1 FIRM / 2 SOFT / 1 FLOPPY, REDESIGN | SHAPE ρ 0.29–0.55 (4 of 5 sig.); genome control 43rd–82nd %ile, p 0.18–0.57 | FIRM retention 0.50 Andronescu/Langdon, 1.00 at all temps/flanks; ensemble defect 14.68 (37°C) vs 13.55 (42°C) | ViennaRNA cannot represent the FSE pseudoknot; tiers describe a pseudoknot-free approximation | [layer4_shape.md](docs/benchmark/layer4_shape.md) |
| **SMN2 ISS-N1** | 2 FIRM / 6 SOFT / 1 FLOPPY, REDESIGN | — | FIRM retention 0.00 under Langdon2018, 0.40 with 25–100 nt genomic flanks | No SHAPE data for this locus | [layer5_robustness.md](docs/benchmark/layer5_robustness.md) |
| **CFTR 5'UTR** | 3 FIRM / 1 SOFT / 6 FLOPPY, REDESIGN | — | 3 FIRM stems retained under every temperature, parameter set and flank; window = 5′UTR nt 1–70 + first 130 nt CDS | Mostly FLOPPY; not a pure 5′UTR window | [layer0_cases.md](docs/benchmark/layer0_cases.md) |
| **MAPT exon 10** | 5 FIRM / 2 SOFT / 4 FLOPPY, REDESIGN | — | FIRM retention 0.29 at 25-nt flanks, 1.00 at 50/100 nt | FIRM stems sensitive to 25-nt genomic flanks (retention 0.29); stable at 50/100 nt | [layer5_robustness.md](docs/benchmark/layer5_robustness.md) |
| **HCV IRES domain II** | 2 FIRM / 3 SOFT, NEED PROBING | — | Every FIRM and SOFT stem lost with any genomic flank (retention 0.00 at 25/50/100 nt) | Short window; the predicted domain II stems are not stable once genomic flanks are added (retention 0.00), although domain II itself is an NMR-determined structure (Lukavsky 2003) | [layer5_robustness.md](docs/benchmark/layer5_robustness.md) |

**Note:** Tier summary column = case_metrics.csv stem tiers (minimum pair probability per stem); Layer 5 counts use mean-probability stem tiers (layer5 CSVs), so stem counts differ (e.g. FSE 2 FIRM / 3 SOFT, HCV 3 / 3, MAPT 7 FIRM / 1 SOFT).

---

## Limitations

1. **Thermodynamic model only:** ViennaRNA predicts MFE/ensemble structures; does not incorporate phylogenetic covariation or crystallographic constraints. Disagreement with comparative structures is expected and does not invalidate tier reliability for *predicted* stems.

2. **Length caps:** ≤500 nt for ArchiveII/bpRNA/Rfam (Rfam seed projection additionally excluded >400 nt). Full-length 16S/23S rRNA not included.

3. **SHAPE coverage:** Layer 4 tests only the SARS-CoV-2 FSE. Other disease cases lack experimental probing data.

4. **FLOPPY tier uninformative:** PPV ≈ 0.15 means FLOPPY pairs rarely match references. This is correct behavior (low probability = low reliability), not a bug. Users should not trust FLOPPY predictions.

5. **Parameter sensitivity:** Andronescu2007/Langdon2018 cause substantial tier shifts (Layer 5). Results are specific to Turner2004 (ViennaRNA default).

6. **Context dependence:** Genomic flanks alter predictions (Layer 5). Disease windows are isolated sequences; in vivo context may differ.

7. **No pseudoknots:** ViennaRNA does not predict pseudoknots; removed from references before comparison.

---

## Reproduce

### Fetch and verify data bundle

```bash
# Download data bundle (instructions in MANIFEST.md)
# Or extract from attached uploads:
mkdir -p data/_cache
tar -xzf uploads/foldtrust_bench_data.tar.gz -C data/_cache
tar -xzf uploads/bprna_TS0_canonicals.tar.gz -C data/_cache

# Verify bundle
python3 scripts/verify_bundle.py
```

### Run all layers

```bash
# Default subsets (~3 minutes on 16GB machine)
foldtrust benchmark all --output benchmarks/outputs

# Full datasets (longer, adds --full flag)
foldtrust benchmark all --output benchmarks/outputs_full --full
```

### Individual layers

```bash
foldtrust benchmark layer1 --output benchmarks/outputs
foldtrust benchmark layer2 --output benchmarks/outputs
foldtrust benchmark layer3 --output benchmarks/outputs
foldtrust benchmark layer4_shape --output benchmarks/outputs
foldtrust benchmark layer5 --output benchmarks/outputs
```

All outputs are saved to CSV/JSON in `benchmarks/outputs/layerN/`. Figures in `benchmarks/outputs/figures/`.

---

## References

**Data sources** (see `MANIFEST.md` and `SHA256SUMS` for full citations and file hashes):
- ArchiveII: Sloma & Mathews, RNA 2016
- Rfam 15: Ontiveros-Palacios et al., NAR 2025 (release 15.1)
- bpRNA-1m: Danaee et al., NAR 2018; Singh et al., Nat Commun 2019
- SARS-CoV-2 SHAPE: DasLab/SARS_CoV-2_shape_comparison (Incarnato, Pyle, Zhang labs)
- ViennaRNA 2.7.x: Lorenz et al., Algorithms Mol Biol 2011

**Detailed layer documentation:** See `docs/benchmark/layer*.md` for methods, results, and interpretation.

**Code:** All benchmark code in `src/foldtrust/benchmark/`, tests in `tests/`.
