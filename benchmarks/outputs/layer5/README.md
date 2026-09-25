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

**Runtime:** ~10-15 seconds on a 16GB Mac (5 cases × 3 temperatures × 2 param sets × 4 window contexts = 100+ folds)

## Data Sources

All five cases use Layer 0 verified coordinates:
- `sars2-fse`: NC_045512.2:13462-13542 (81 nt)
- `smn2-iss-n1`: NG_008728.1:31999-32152 (154 nt)
- `cftr-5utr`: NM_000492.4:1-200 (200 nt)
- `mapt-e10`: NG_007398.2:120818-121000 (183 nt)
- `hcv-ires-dii`: AF009606.1:44-118 (75 nt)

Flanking sequences fetched from NCBI and cached in `data/_cache/` (not committed to repo).

## Key Findings

1. **FIRM stems** (p ≥ 0.85): 100% retained across temperatures, 80-88% across parameter sets, 65% with 100 nt flanks
2. **SOFT stems** (0.5 ≤ p < 0.85): 60-65% retained across temperatures, 44-52% across parameter sets, 44% with 100 nt flanks
3. **Parameter sets matter:** Andronescu2007 changes MFE by 3-5 kcal/mol per case vs Turner2004
4. **Window context is critical:** Adding 100 nt flanks increases BP distance by ~80-90 pairs

## Tests

See `tests/test_benchmark_layer5.py`:
- Case coordinate verification (match Layer 0)
- FSE 37°C Turner2004 energy = -26.00 kcal/mol
- Self-retention = 1.0
- Flank extraction verified as exact substring

All tests pass: `pytest tests/test_benchmark_layer5.py -v`
