# Layer 5: Robustness Analysis Outputs

This directory contains all results from the Layer 5 robustness benchmark.

## Files

### Detailed Results (JSON)
- **`layer5_complete.json`** — Complete results for all cases, all conditions (temperature, parameters, window context)

### Temperature Sweep (Turner2004 at 25, 30, 42°C vs 37°C baseline)
- **`temperature_mfe_energies.csv`** — MFE energies at each temperature
- **`temperature_bp_distances.csv`** — Base-pair distances (MFE and MEA to baseline), ensemble defects
- **`temperature_stem_retention.csv`** — Stem retention by tier (FIRM/SOFT/FLOPPY) and mean probability change

### Parameter Set Sweep (Andronescu2007, Langdon2018 vs Turner2004 at 37°C)
- **`parameters_mfe_energies.csv`** — MFE energies for Turner2004, Andronescu2007, Langdon2018
- **`parameters_bp_distances.csv`** — Base-pair distances and ensemble defects
- **`parameters_stem_retention.csv`** — Stem retention by tier

### Window Context Sweep (0, 25, 50, 100 nt flanks, Turner2004 at 37°C)
- **`window_context_metrics.csv`** — Extended sequence lengths, MFE energies, BP distances
- **`window_context_stem_retention.csv`** — Core-window stem retention when flanking sequence added

### Summary Tables
- **`layer5_tables.md`** — Rendered summary tables (pooled and per-case mean retention) from CSVs

### Figures
Generated figures are in `benchmarks/outputs/figures/`:
- `layer5_temperature_retention.png` — Stem retention heatmap across temperatures
- `layer5_params_retention.png` — Stem retention heatmap across parameter sets
- `layer5_window_retention.png` — Core-stem retention vs flank size

### Decision Document
- **`VIENNA_DECISION.md`** — Rationale for subprocess-based parameter isolation in `vienna.py`

## Command to Reproduce

```python
from pathlib import Path
from foldtrust.benchmark.layer5_robustness import run_layer5_analysis

results = run_layer5_analysis(
    cases_dir=Path('data/cases'),
    output_dir=Path('benchmarks/outputs/layer5'),
    cache_dir=Path('data/_cache')
)
```

Render summary tables:
```bash
python3 scripts/render_layer5_tables.py
```

**Runtime:** ~10-15 seconds (5 cases × 3 temperatures × 2 param sets × 4 window contexts = 100+ folds)

## Data Sources

All five cases use Layer 0 verified coordinates:
- `sars2-fse`: NC_045512.2:13462-13542 (81 nt)
- `smn2-iss-n1`: NG_008728.1:31999-32152 (154 nt)
- `cftr-5utr`: NM_000492.4:1-200 (200 nt)
- `mapt-e10`: NG_007398.2:120818-121000 (183 nt)
- `hcv-ires-dii`: AF009606.1:44-118 (75 nt)

Flanking sequences fetched from NCBI and cached in `data/_cache/` (not committed to repo).

**Flank clipping:**
- CFTR left flank: 0 nt (clipped at record start)
- HCV left flank: 43 nt at 50/100 nt flank sizes (clipped at record start, position 1)
- All other cases: no clipping

## Key Findings (Pooled Retention)

All numbers from `layer5_tables.md` (computed from CSVs).

### Baseline (37°C Turner2004 MFE stems vs MEA structure)
- **FIRM** (p ≥ 0.85): 1.000 (100% of MFE FIRM stems retained in MEA)
- **SOFT** (0.5 ≤ p < 0.85): 0.733 (73.3% retained)
- **FLOPPY** (p < 0.5): 0.364 (36.4% retained)

The baseline SOFT retention of 0.733 reflects the **MFE-vs-MEA gap**, not a perturbation. Much of the SOFT "loss" at other conditions is this same gap.

### Temperature (25, 30, 42°C vs 37°C baseline)
- **FIRM**: 95% retained at 25/30°C, 100% at 42°C
- **SOFT**: 60% retained at 25/30°C (pooled), 73.3% at 42°C (matches baseline)
- **FLOPPY**: 27-36% (baseline already low)

### Parameters (Andronescu2007, Langdon2018 vs Turner2004)
- **FIRM**: 80% under Andronescu2007, 60% under Langdon2018
- **SOFT**: 6.7% under Andronescu2007, 13.3% under Langdon2018
- **MFE energy shifts**: Andronescu2007 +2.23 to +4.91 kcal/mol vs Turner2004 (per case), Langdon2018 -1.40 to +2.20 kcal/mol

### Window Context (25, 50, 100 nt flanks vs 0 nt)
- **FIRM**: 45% with 25 nt flanks, 70% with 50/100 nt flanks
- **SOFT**: 20% with 25/50 nt flanks, 13.3% with 100 nt flanks
- **Mean MEA BP distance** (core pairs only): 0 nt → 30.4 bp (25 nt), 19.6 bp (50 nt), 22.4 bp (100 nt)

### Ensemble Defect (42°C vs 37°C)
Only SMN2 shows increased defect at 42°C (+9.5%). Four of five cases show decreased or stable defect (-23.7% to -2.5%). The claim "42°C ensemble defect 10-15% higher" is false.

## Method Notes

- **MFE stems** parsed by `vienna.parse_stems`: maximal stacked runs ≥2 pairs, tiered by mean stem pair probability (FIRM ≥0.85, SOFT 0.5-0.85, FLOPPY <0.5)
- **Condition structure** = ViennaRNA `fc.MEA(gamma=1)` after `fc.pf()`
- **Stem is retained** if all its pairs are present in the condition MEA
- **Flank BP distance** counts only pairs with both ends in the core (neither end in flanking sequence)
- **Parameter sets** run in separate subprocesses (ViennaRNA 2.7.2 `params_load` + fresh `md` does not take effect in a process that already folded)

## Parameter Set Citations

- **Andronescu2007**: Andronescu, Condon, Hoos, Mathews & Murphy 2007, "Efficient parameter estimation for RNA secondary structure prediction", *Bioinformatics* 23:i19-i28, DOI [10.1093/bioinformatics/btm223](https://doi.org/10.1093/bioinformatics/btm223)
- **Langdon2018**: Langdon, Petke & Lorenz 2018, "Evolving better RNAfold structure prediction" (EuroGP 2018, genetic improvement of RNAfold energy parameters), DOI [10.1007/978-3-319-77553-1_14](https://doi.org/10.1007/978-3-319-77553-1_14)

## Tests

See `tests/test_benchmark_layer5.py`:
- Self-retention = 1.0 (MFE stems retained in same MFE structure)
- Baseline retention = saved baseline row (37°C Turner2004 MFE stems in 37°C Turner2004 MEA)
- Flank extraction: left flank ends exactly where core starts, right flank begins exactly where core ends, at recorded coordinates, in cached NCBI record

All tests pass: `pytest tests/test_benchmark_layer5.py -v`
