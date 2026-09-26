# Layer 5: Robustness — Structure Stability Across Conditions

## Question

How stable are FoldTrust's MFE structure, MEA structure, and tier classifications when biological or methodological conditions vary?

RNA structure prediction relies on thermodynamic models (parameter sets), temperature assumptions, and window boundary choices. Layer 5 tests whether FoldTrust's predictions and reliability tiers remain consistent under:

1. **Temperature variation** — 25, 30, 37, 42°C
2. **Energy parameter sets** — Turner2004 (default), Andronescu2007, Langdon2018
3. **Window context** — extending the sequence window with 0, 25, 50, 100 nt of real flanking sequence from the same NCBI record

All conditions are compared to the **baseline**: Turner2004 parameters, 37°C, core window only (0 nt flanks).

## Data

The five Layer 0 disease-relevant cases (verified coordinates):

| Case | Accession | Coordinates | Length | Description |
|------|-----------|-------------|--------|-------------|
| **sars2-fse** | NC_045512.2 | 13462-13542 | 81 nt | SARS-CoV-2 frameshift element |
| **smn2-iss-n1** | NG_008728.1 | 31999-32152 | 154 nt | SMN2 ISS-N1 |
| **cftr-5utr** | NM_000492.4 | 1-200 | 200 nt | CFTR 5' UTR + start codon |
| **mapt-e10** | NG_007398.2 | 120818-121000 | 183 nt | MAPT exon 10 splice site |
| **hcv-ires-dii** | AF009606.1 | 44-118 | 75 nt | HCV IRES domain II |

## Method

### Baseline (Reference Condition)

For each case:
1. Fold with **Turner2004** parameters at **37°C** (core window only)
2. Compute MFE structure using ViennaRNA
3. Compute MEA structure using `fc.MEA(gamma=1.0)` after `fc.pf()`
4. Parse stems from MFE structure using `vienna.parse_stems`:
   - Stem = maximal run of ≥2 stacked pairs (i,j), (i+1,j-1), ...
   - Stem tier based on **mean pair probability** of stem's pairs:
     - **FIRM**: mean p ≥ 0.85
     - **SOFT**: 0.5 ≤ mean p < 0.85
     - **FLOPPY**: mean p < 0.5

### Perturbation Tests

#### 1. Temperature Sweep
Fold each case at **25, 30, 42°C** with Turner2004, compare to 37°C baseline.

**Metrics per condition:**
- MFE energy (kcal/mol)
- Ensemble defect (expected number of incorrectly paired nucleotides)
- Base-pair distance of MFE structure to baseline 37°C MEA
- Base-pair distance of condition MEA to baseline 37°C MEA
- Stem retention by tier (fraction of baseline stems whose pairs are all paired in condition MEA)
- Mean change in stem mean probability

#### 2. Parameter Set Sweep
Fold each case with **Andronescu2007** and **Langdon2018** parameters at 37°C, compare to Turner2004.

**Metrics:** Same as temperature sweep.

**Parameter Set Details:**
- **Turner2004**: Standard ViennaRNA parameters (nearest-neighbor model)
- **Andronescu2007**: Andronescu, Condon, Hoos, Mathews & Murphy 2007, "Efficient parameter estimation for RNA secondary structure prediction", *Bioinformatics* 23:i19-i28, DOI [10.1093/bioinformatics/btm223](https://doi.org/10.1093/bioinformatics/btm223)
- **Langdon2018**: Langdon, Petke & Lorenz 2018, "Evolving better RNAfold structure prediction" (EuroGP 2018, genetic improvement of RNAfold energy parameters), DOI [10.1007/978-3-319-77553-1_14](https://doi.org/10.1007/978-3-319-77553-1_14)

ViennaRNA parameter state is isolated via subprocess calls because `params_load()` + fresh `md` does not take effect in a process that already folded. See `benchmarks/outputs/layer5/VIENNA_DECISION.md`.

#### 3. Window Context Sweep
Extend each case window by adding **0, 25, 50, 100 nt** of real flanking sequence from the same NCBI accession, fold the extended sequence with Turner2004 at 37°C, and test **retention of core-window stems**.

**Flank extraction:**
- Fetch complete record from NCBI (cached in `data/_cache/`)
- Extract left flank: positions `[start - flank_size, start - 1]`
- Extract right flank: positions `[end + 1, end + flank_size]`
- Clip at record boundaries
- Verify core sequence is exact substring of cached record

**Flank coordinates and clipping** (from `window_context_metrics.csv`):

| Case | Flank Size | Left Coords | Right Coords | Left Clipped | Right Clipped |
|------|------------|-------------|--------------|--------------|---------------|
| sars2-fse | 25 | NC_045512.2:13437-13461 | NC_045512.2:13543-13567 | No | No |
| sars2-fse | 50 | NC_045512.2:13412-13461 | NC_045512.2:13543-13592 | No | No |
| sars2-fse | 100 | NC_045512.2:13362-13461 | NC_045512.2:13543-13642 | No | No |
| smn2-iss-n1 | 25 | NG_008728.1:31974-31998 | NG_008728.1:32153-32177 | No | No |
| smn2-iss-n1 | 50 | NG_008728.1:31949-31998 | NG_008728.1:32153-32202 | No | No |
| smn2-iss-n1 | 100 | NG_008728.1:31899-31998 | NG_008728.1:32153-32252 | No | No |
| cftr-5utr | 25 | (none) | NM_000492.4:201-225 | Yes (0 nt) | No |
| cftr-5utr | 50 | (none) | NM_000492.4:201-250 | Yes (0 nt) | No |
| cftr-5utr | 100 | (none) | NM_000492.4:201-300 | Yes (0 nt) | No |
| mapt-e10 | 25 | NG_007398.2:120793-120817 | NG_007398.2:121001-121025 | No | No |
| mapt-e10 | 50 | NG_007398.2:120768-120817 | NG_007398.2:121001-121050 | No | No |
| mapt-e10 | 100 | NG_007398.2:120718-120817 | NG_007398.2:121001-121100 | No | No |
| hcv-ires-dii | 25 | AF009606.1:19-43 | AF009606.1:119-143 | No | No |
| hcv-ires-dii | 50 | AF009606.1:1-43 | AF009606.1:119-168 | Yes (43 nt) | No |
| hcv-ires-dii | 100 | AF009606.1:1-43 | AF009606.1:119-218 | Yes (43 nt) | No |

**Retention metric:**
For each baseline stem (FIRM/SOFT/FLOPPY), count whether **all** its pairs are paired in the extended-window MEA. A stem is retained if all its pairs are present in the condition MEA structure. Only core-window stems are counted; new stems in flanking regions are ignored.

### Ensemble Defect Calculation

Ensemble defect = Σ (1 - p_correct) over all positions, where:
- If position i is paired to j in the structure: p_correct = p(i,j)
- If position i is unpaired: p_correct = 1 - Σ_k p(i,k)

### Base-Pair Distance

|structure1 - structure2| = number of pairs in symmetric difference = |(pairs1 \ pairs2)| + |(pairs2 \ pairs1)|

### Flank BP Distance

For window context tests, `bp_distance_mea` counts only pairs with both ends in the core window (neither end in the flanking sequence). This isolates the structural change within the region of interest.

## Command

```bash
# Run complete Layer 5 analysis
python3 -c "
from pathlib import Path
from foldtrust.benchmark.layer5_robustness import run_layer5_analysis

results = run_layer5_analysis(
    cases_dir=Path('data/cases'),
    output_dir=Path('benchmarks/outputs/layer5'),
    cache_dir=Path('data/_cache')
)
"
```

Output CSVs:
- `temperature_mfe_energies.csv` — MFE energies at each temperature
- `temperature_bp_distances.csv` — BP distances and ensemble defects
- `temperature_stem_retention.csv` — Stem retention by tier and temperature
- `parameters_mfe_energies.csv` — MFE energies for each parameter set
- `parameters_bp_distances.csv` — BP distances and ensemble defects by parameter set
- `parameters_stem_retention.csv` — Stem retention by tier and parameter set
- `window_context_metrics.csv` — Extended sequence lengths, MFE energies, BP distances
- `window_context_stem_retention.csv` — Core-stem retention with flanking context
- `layer5_complete.json` — Full detailed results

Render summary tables:
```bash
python3 scripts/render_layer5_tables.py
```

Output: `benchmarks/outputs/layer5/layer5_tables.md`

## Results

All numbers below are from the CSVs and rendered tables.

### MFE Energies by Parameter Set

| case | Turner2004_baseline | Andronescu2007 | Langdon2018 |
|------|---------------------|----------------|-------------|
| sars2-fse | -26.00 | -22.26 | -24.70 |
| smn2-iss-n1 | -31.30 | -26.39 | -29.10 |
| cftr-5utr | -60.60 | -56.03 | -62.00 |
| mapt-e10 | -52.30 | -50.07 | -51.80 |
| hcv-ires-dii | -23.80 | -21.26 | -22.00 |

**Energy shifts:** Andronescu2007 shifts MFE by +2.23 to +4.91 kcal/mol (less stable) vs Turner2004. Langdon2018 shifts by -1.40 to +2.20 kcal/mol.

**SARS-CoV-2 FSE verification:** Turner2004 at 37°C gives -26.00 kcal/mol (matches Layer 0 and `tests/test_param_energies.py`).

### Temperature Sweep: Stem Retention

Pooled retention (total retained stems / total reference stems across all 5 cases):

| tier   | 37C_baseline | 25.0C | 30.0C | 42.0C |
|:-------|-------------:|------:|------:|------:|
| FIRM   | 1.000 | 0.950 | 0.950 | 1.000 |
| SOFT   | 0.733 | 0.600 | 0.600 | 0.733 |
| FLOPPY | 0.364 | 0.273 | 0.364 | 0.364 |

**Change vs baseline:**
- **FIRM** (p ≥ 0.85): 25/30°C lose 5% (pooled); 42°C retains 100%
- **SOFT** (0.5 ≤ p < 0.85): 25/30°C retain 60% (pooled), 42°C matches baseline
- **FLOPPY** (p < 0.5): baseline retention already low (36.4%), minimal change

**Note:** The baseline SOFT retention of 0.733 (not 1.0) reflects the **MFE-vs-MEA gap**: 73.3% of Turner2004 37°C MFE SOFT stems are retained in the Turner2004 37°C MEA. Much of the apparent "loss" at other temperatures is this same MFE-MEA gap, not the temperature perturbation itself. The baseline row is computed by comparing MFE stems (parsed from MFE structure) to MEA stems (parsed from MEA structure) at the same condition.

### Parameter Set Sweep: Stem Retention

Pooled retention (total retained stems / total reference stems):

| tier   | Turner2004_baseline | Andronescu2007 | Langdon2018 |
|:-------|--------------------:|---------------:|------------:|
| FIRM   | 1.000 | 0.800 | 0.600 |
| SOFT   | 0.733 | 0.067 | 0.133 |
| FLOPPY | 0.364 | 0.455 | 0.091 |

**Change vs baseline:**
- **FIRM**: Andronescu2007 retains 80% (-20%), Langdon2018 retains 60% (-40%)
- **SOFT**: Andronescu2007 retains 6.7% (-66.6 pp), Langdon2018 retains 13.3% (-60.0 pp)
- **FLOPPY**: Andronescu2007 retains 45.5% (+9.1 pp), Langdon2018 retains 9.1% (-27.3 pp)

**Observation:** Parameter set changes affect even high-confidence stems. SOFT stems are especially sensitive.

### Window Context: Core-Stem Retention

Pooled retention of core-window stems when flanking sequence is added:

| tier   | 0 nt | 25 nt | 50 nt | 100 nt |
|:-------|-----:|------:|------:|-------:|
| FIRM   | 1.000 | 0.450 | 0.700 | 0.700 |
| SOFT   | 0.733 | 0.200 | 0.200 | 0.133 |
| FLOPPY | 0.364 | 0.545 | 0.364 | 0.455 |

Mean MEA base-pair distance (core pairs only):

| Flank Size (nt) | Mean MEA BP Distance |
|-----------------|----------------------|
| 0 | 0.0 |
| 25 | 30.4 |
| 50 | 19.6 |
| 100 | 22.4 |

**Change vs baseline:**
- **FIRM**: 25 nt flanks drop retention to 45.0% (-55.0 pp), 50/100 nt to 70.0% (-30.0 pp)
- **SOFT**: 25/50 nt flanks retain 20.0% (-53.3 pp), 100 nt retains 13.3% (-60.0 pp)
- **FLOPPY**: 25 nt flanks increase retention to 54.5% (+18.2 pp)

**Observation:** Adding flanking sequence dramatically reduces FIRM and SOFT stem retention, because competing structures in the extended sequence can sequester nucleotides. FLOPPY stems sometimes increase retention (they were already unstable in the baseline).

### Ensemble Defect: Temperature Effects

Ensemble defect change at 42°C vs 37°C baseline:

| Case | 37°C Defect | 42°C Defect | Change | % Change |
|------|-------------|-------------|--------|----------|
| sars2-fse | 14.68 | 13.55 | -1.13 | -7.7% |
| smn2-iss-n1 | 22.26 | 24.37 | +2.11 | +9.5% |
| cftr-5utr | 54.98 | 52.18 | -2.80 | -5.1% |
| mapt-e10 | 32.99 | 25.16 | -7.83 | -23.7% |
| hcv-ires-dii | 10.92 | 10.65 | -0.27 | -2.5% |

**Observation:** Only SMN2 shows increased ensemble defect at 42°C (+9.5%). Four of five cases show decreased or stable defect.

### Figures

Generated figures in `benchmarks/outputs/figures/`:
- `layer5_temperature_retention.png` — Stem retention heatmap across temperatures
- `layer5_params_retention.png` — Stem retention heatmap across parameter sets
- `layer5_window_retention.png` — Core-stem retention vs flank size

## Interpretation

1. **FIRM stems are temperature-robust but parameter-sensitive.** 95% of FIRM stems (pooled) are retained from 25-30°C, 100% at 42°C. But only 80% are retained under Andronescu2007, and 60% under Langdon2018.

2. **SOFT stems show moderate stability.** 60% are retained at 25/30°C (pooled). The baseline SOFT retention of 73.3% reflects the MFE-vs-MEA gap: even at the same condition, not all MFE SOFT stems are retained in the MEA structure. This gap is present across all conditions.

3. **Window context strongly affects predictions.** Even FIRM stems drop to 45-70% retention when 25-100 nt flanks are added, because competing structures in the extended sequence can sequester nucleotides.

4. **Ensemble defect does not consistently increase with temperature.** At 42°C, only SMN2 shows increased defect (+9.5%); the other four cases show decreased or stable defect (-23.7% to -2.5%).

5. **Base-pair distance correlates with window extension.** Adding 25-100 nt flanks changes ~20-30 base pairs (core pairs only) on average across the five cases.

## Limitations

1. **MEA structure may not represent biological ground truth.** ViennaRNA's MEA maximizes expected accuracy but does not guarantee the most abundant structure in vivo. Boltzmann sampling or kinetic folding may reveal additional structural heterogeneity.

2. **Temperature range tested is narrow.** 25-42°C covers room temperature to mild fever, but does not test extreme denaturation or cold-induced refolding.

3. **Parameter sets tested are MFE-focused.** Other models (e.g., CONTRAfold, RNAstructure's partition function parameters) were not included.

4. **Flanking sequences assume no trans interactions.** Real cellular context includes RNA-binding proteins, ribosomes, and other RNAs that may stabilize or disrupt local structure.

5. **Retention metric is binary.** A stem is either fully retained or not; partial retention (e.g., 3 of 4 pairs) is counted as not retained.

6. **No experimental validation.** SHAPE, DMS-MaPseq, or NMR data at different temperatures or with extended windows would be needed to validate these predictions.

## Citations

- Andronescu M, Condon A, Hoos HH, Mathews DH, Murphy KP. 2007. Efficient parameter estimation for RNA secondary structure prediction. *Bioinformatics* 23(13):i19-i28. DOI: [10.1093/bioinformatics/btm223](https://doi.org/10.1093/bioinformatics/btm223)
- Langdon WB, Petke J, Lorenz R. 2018. Evolving better RNAfold structure prediction. In: *Lecture Notes in Computer Science* vol 10781, pp 220-236. Springer. DOI: [10.1007/978-3-319-77553-1_14](https://doi.org/10.1007/978-3-319-77553-1_14)
- Lorenz R, Bernhart SH, Höner zu Siederdissen C, Tafer H, Flamm C, Stadler PF, Hofacker IL. 2011. ViennaRNA Package 2.0. *Algorithms Mol Biol* 6:26. DOI: [10.1186/1748-7188-6-26](https://doi.org/10.1186/1748-7188-6-26)
