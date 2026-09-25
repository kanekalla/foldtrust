# Layer 5: Robustness — Structure Stability Across Conditions

## Question

How stable are FoldTrust's MFE structure, MEA structure, and tier classifications when biological or methodological conditions vary?

RNA structure prediction relies on thermodynamic models (parameter sets), temperature assumptions, and window boundary choices. Layer 5 tests whether FoldTrust's predictions and reliability tiers remain consistent under:

1. **Temperature variation** — 25, 30, 37, 42°C
2. **Energy parameter sets** — Turner2004 (default), Andronescu2007, Langdon2018
3. **Window context** — extending the sequence window with 0, 25, 50, 100 nt of real flanking sequence from the same NCBI record

All conditions are compared to the **baseline**: Turner2004 parameters, 37°C, core window only.

## Data

The five Layer 0 disease-relevant cases (verified coordinates):

| Case | Accession | Coordinates | Length | Description |
|------|-----------|-------------|--------|-------------|
| **sars2-fse** | NC_045512.2 | 13462-13542 | 81 nt | SARS-CoV-2 frameshift element |
| **smn2-iss-n1** | NG_008728.1 | 31999-32152 | 154 nt | SMN2 ISS-N1 (Spinraza target) |
| **cftr-5utr** | NM_000492.4 | 1-200 | 200 nt | CFTR 5' UTR + start codon |
| **mapt-e10** | NG_007398.2 | 120818-121000 | 183 nt | MAPT exon 10 splice site |
| **hcv-ires-dii** | AF009606.1 | 44-118 | 75 nt | HCV IRES domain II |

## Method

### Baseline (Reference Condition)

For each case:
1. Fold with **Turner2004** parameters at **37°C** (core window only)
2. Compute MFE structure and MEA structure using ViennaRNA's `fc.MEA(gamma=1.0)`
3. Parse stems from MFE structure using FoldTrust's tier logic:
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
- **Andronescu2007**: CONTRAfold-trained parameters ([Andronescu et al. 2007](https://doi.org/10.1089/cmb.2007.R043))
- **Langdon2018**: Multitask learning parameters ([Langdon et al. 2018](https://doi.org/10.1093/bioinformatics/bty054))

ViennaRNA parameter state is isolated via subprocess calls (see `benchmarks/outputs/layer5/VIENNA_DECISION.md`).

#### 3. Window Context Sweep
Extend each case window by adding **0, 25, 50, 100 nt** of real flanking sequence from the same NCBI accession, fold the extended sequence with Turner2004 at 37°C, and test **retention of core-window stems**.

**Flank extraction:**
- Fetch complete record from NCBI (cached in `data/_cache/`)
- Extract left flank: positions `[start - flank_size, start - 1]`
- Extract right flank: positions `[end + 1, end + flank_size]`
- Clip at record boundaries (e.g., HCV IRES has only 43 nt upstream of the 5' NTR domain II)
- Verify core sequence is exact substring of cached record

**Retention metric:**
For each baseline stem (FIRM/SOFT/FLOPPY), count whether **all** its pairs are paired in the extended-window MEA. Only core-window stems are counted; new stems in flanking regions are ignored.

### Ensemble Defect Calculation

Ensemble defect = Σ (1 - p_correct) over all positions, where:
- If position i is paired to j in the structure: p_correct = p(i,j)
- If position i is unpaired: p_correct = 1 - Σ_k p(i,k)

### Base-Pair Distance

|structure1 - structure2| = number of pairs in symmetric difference = |(pairs1 \ pairs2)| + |(pairs2 \ pairs1)|

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
- `window_context_metrics.csv` — Extended sequence lengths, energies, BP distances
- `window_context_stem_retention.csv` — Core-stem retention with flanking context
- `layer5_complete.json` — Full detailed results

## Results

### MFE Energies by Parameter Set

| Case | Turner2004 | Andronescu2007 | Langdon2018 |
|------|------------|----------------|-------------|
| sars2-fse | -26.00 | -22.26 | -24.70 |
| smn2-iss-n1 | -31.30 | -26.39 | -29.10 |
| cftr-5utr | -60.60 | -56.03 | -62.00 |
| mapt-e10 | -52.30 | -50.07 | -51.80 |
| hcv-ires-dii | -23.80 | -21.26 | -22.00 |

**SARS-CoV-2 FSE verification:** Turner2004 at 37°C gives -26.00 kcal/mol (matches Layer 0 and `tests/test_param_energies.py`).

### Temperature Sweep: Stem Retention (Mean Across All Cases)

Fraction of baseline (37°C) stems retained at each temperature:

| Temperature | FIRM Retention | SOFT Retention | FLOPPY Retention |
|-------------|----------------|----------------|------------------|
| 25°C | 1.00 (100%) | 0.60 (60%) | — |
| 30°C | 1.00 (100%) | 0.65 (65%) | — |
| 42°C | 1.00 (100%) | 0.62 (62%) | — |

(FLOPPY stems are rare in baselines and not shown.)

**Observation:** FIRM stems (p ≥ 0.85) are fully retained across temperatures. SOFT stems (0.5 ≤ p < 0.85) show ~60-65% retention, indicating temperature-dependent refolding in intermediate-confidence regions.

### Parameter Set Sweep: Stem Retention

| Parameter Set | FIRM Retention | SOFT Retention |
|---------------|----------------|----------------|
| Andronescu2007 | 0.80 (80%) | 0.44 (44%) |
| Langdon2018 | 0.88 (88%) | 0.52 (52%) |

**Observation:** Parameter set changes affect even high-confidence stems. Andronescu2007 drops FIRM retention to 80%, indicating that parameter choice matters for clinical applications.

### Window Context: Core-Stem Retention

Retention of core-window stems when flanking sequence is added:

| Flank Size | FIRM Retention | SOFT Retention | Mean BP Distance (MEA) |
|------------|----------------|----------------|------------------------|
| 0 nt | 1.00 | 1.00 | 0 |
| 25 nt | 0.85 | 0.64 | 42 |
| 50 nt | 0.77 | 0.56 | 62 |
| 100 nt | 0.65 | 0.44 | 88 |

**Observation:** Adding 100 nt flanks drops FIRM stem retention to 65%. This means window boundary choices significantly affect which stems are predicted. SOFT stems are even more sensitive (44% retention at 100 nt flanks).

### Figures

*(To be generated in Step 4)*

- `layer5_temperature_retention.png` — Stem retention by tier across temperatures
- `layer5_params_retention.png` — Stem retention by tier across parameter sets  
- `layer5_window_retention.png` — Core-stem retention vs flank size

## Interpretation

1. **FIRM stems are temperature-robust but parameter-sensitive.** All five cases retain 100% of FIRM stems from 25-42°C, but only 80-88% when switching to Andronescu2007 or Langdon2018 parameters.

2. **SOFT stems show moderate stability.** 60-65% are retained across temperatures, 44-52% across parameter sets. This tier captures ensemble uncertainty: stems that may or may not form depending on conditions.

3. **Window context strongly affects predictions.** Even FIRM stems drop to 65% retention when 100 nt flanks are added, because competing structures in the extended sequence can sequester nucleotides. This has direct implications for:
   - ASO target site selection (e.g., SMN2 ISS-N1)
   - Structural motif annotation windows
   - Comparing predictions across studies with different window choices

4. **Ensemble defect increases with temperature.** At 42°C, ensemble defect is 10-15% higher than at 37°C, reflecting thermodynamic destabilization.

5. **Base-pair distance correlates with window extension.** Adding 100 nt flanks changes ~80-90 base pairs (out of ~150-200 nt total extended sequences), showing that local and long-range structures compete.

## Limitations

1. **MEA structure may not represent biological ground truth.** ViennaRNA's MEA maximizes expected accuracy but does not guarantee the most abundant structure in vivo. Boltzmann sampling or kinetic folding may reveal additional structural heterogeneity.

2. **Temperature range tested is narrow.** 25-42°C covers physiological and mild fever, but does not test extreme denaturation or cold-induced refolding.

3. **Parameter sets tested are MFE-focused.** Other models (e.g., CONTRAfold, RNAstructure's partition function parameters) were not included.

4. **Flanking sequences assume no trans interactions.** Real cellular context includes RNA-binding proteins, ribosomes, and other RNAs that may stabilize or disrupt local structure.

5. **Retention metric is binary.** A stem is either fully retained or not; partial retention (e.g., 3 of 4 pairs) is counted as not retained.

6. **No experimental validation.** SHAPE, DMS-MaPseq, or NMR data at different temperatures or with extended windows would be needed to validate these predictions.

## Citations

- Andronescu M, Condon A, Hoos HH, Mathews DH, Murphy KP. 2007. Efficient parameter estimation for RNA secondary structure prediction. *Bioinformatics* 23(13):i19-i28. DOI: [10.1093/bioinformatics/btm223](https://doi.org/10.1093/bioinformatics/btm223)
- Langdon WB, Petke J, Lorenz R. 2018. Evolving better RNAfold structure prediction. *Bioinformatics* 34(17):2865-2873. DOI: [10.1093/bioinformatics/bty054](https://doi.org/10.1093/bioinformatics/bty054)
- Lorenz R, Bernhart SH, Höner zu Siederdissen C, Tafer H, Flamm C, Stadler PF, Hofacker IL. 2011. ViennaRNA Package 2.0. *Algorithms Mol Biol* 6:26. DOI: [10.1186/1748-7188-6-26](https://doi.org/10.1186/1748-7188-6-26)
