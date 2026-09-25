# Layer 1-3 Fixes Progress

## Completed Fixes

### Fix #2: Restore fetch_benchmark_data.py ✓
- Commit 108912d: Restored from commit 3b89ba6
- Mode set to 100755
- Full fetch script with SHA256 verification

### Fix #3: Remove tRNA-Phe ✓
- Commit 48309bf: Removed tRNA_phe_yeast.{ct,fa,txt}
- Removed entry from curated_references.json
- Fixed docs/API.md and src/foldtrust/io.py examples to use 5S_rRNA_ecoli
- Reason: Reference was ViennaRNA MFE prediction, 75-char dot-bracket for 76-nt sequence

### Fix #6: Slip-tolerant scoring ✓
- Commit f0384cf: Replaced greedy two-sided slip implementation
- New implementation: Standard one-side ±1
  - PPV_slip = #pred with match / #pred
  - SEN_slip = #ref recovered / #ref
  - Match: (i,j), (i±1,j), (i,j±1)
- Added test asserting slip ≥ exact

### Fix #12: tl.py ensemble bugs ✓  
- Commit de81641 (local, not pushed separately): Fixed compute_ensemble
  - Changed bpp copy from 0-based to correct 1-based ViennaRNA indexing
  - Fill symmetric matrix (both triangles)
  - Fixed compute_unpaired_probs: use 1 - P.sum(axis=1), not divide by 2

### Dependencies ✓
- Added ViennaRNA>=2.7 to pyproject.toml

### File structure ✓
- Created scripts/run_layers_1_2_3.py runner
- Created tests/test_benchmark_layer1.py
- Created tests/test_tl_ensemble_regression.py

## Remaining Fixes (TODO)

### Fix #4: Single sampling function
- Layer 2 currently re-seeds before each dataset (lines 365, 373, 381)
- Need to save sample_ids.csv with dataset,name,length,seq_sha1
- Layer 3 must read that file
- KEEP Layer 2 selection (same IDs as saved layer2_raw_results.csv)

### Fix #5: Seed every bootstrap
- Pass `np.random.default_rng(seed)` to bootstrap functions
- Record seed and B in outputs

### Fix #7: Tier PPV
- Headline must be POOLED PPV per tier over MFE pairs
- Add macro mean and N (structs with ≥1 pair in tier)
- Add MEA pair analysis
- Add all-candidate tier PPV (p>1e-3)
- Add stem-level tiers matching FoldTrust parse_stems logic
- State all definitions in docs

### Fix #8: Calibration wording
- Remove "well-calibrated", "validates the core claim", base-rate "5.9×/1.6×"
- Say plainly: probabilities rank well (AUROC ~0.887) but OVERCONFIDENT for p≥0.5
- Include full per-bin table and ECE restricted to p≥0.5
- Note ECE small only because ~74% in [0,0.1)

### Fix #9: Generate numbers from scripts
- Create scripts/render_layer123_tables.py
- Generate all markdown tables from CSV/JSON
- Fix known wrong numbers (see task lines 47-52)
- Fix citations and DOIs (task lines 53-63)
- Fix benchmarks/outputs/README.md

### Fix #10: State exact matching in layer2 docs
- Headline is EXACT pair matching, F1 = mean of per-structure F1
- Reference PKs removed with greedy rule
- Report pooled F1 and slip F1 per method and dataset

### Fix #11: Layer 1 as real pytest
- Create tests/test_benchmark_layer1.py with all coverage
- Include energy regression, bpp symmetry, unpaired prob, GC hairpin FIRM
- Write per-test results to layer1_tests.json AND CSV

### Fix #13: Wire CLI subcommand
- `foldtrust benchmark layer1|layer2|layer3 [--full]`
- Expose as ft.tl.benchmark.layer1/2/3
- Move plotting to ft.pl, make colourblind-safe
- Regenerate figures

### Final Steps
- Re-run default benchmark: `python scripts/run_layers_1_2_3.py`
- Save log to benchmarks/outputs/layer123_run.log
- Regenerate doc tables with render script
- Run pytest (must pass)
- Run ruff/black on changed files
- Verify expected values match (see task lines 78-90)

## Expected Values (for verification)
Same 600 Layer-2 IDs, Turner2004 37°C, exact matching:
- Overall MFE F1 0.548 (sen 0.639, ppv 0.499)
- MEA F1 0.563 (0.641/0.521)
- Centroid F1 0.570 (0.622/0.551)
- Per dataset: ArchiveII 0.562/0.582/0.591, TS0 0.514/0.525/0.532, Rfam 0.569/0.582/0.586
- Slip F1: MFE 0.579, MEA 0.594, centroid 0.599

Layer 3:
- ECE 0.0664, AUROC 0.8889, AUPRC 0.6153
- ECE (p≥0.5) 0.317
- MFE-pair pooled: FIRM 0.674, SOFT 0.299, FLOPPY 0.146
- See full per-bin table in task lines 83-86

## Git Status
Commits pushed to origin/main:
- 108912d: Restore fetch script
- 48309bf: Remove tRNA-Phe  
- f0384cf: Fix slip scoring and add ViennaRNA

NOTE: Co-authored-by trailers present in pushed commits (could not strip them after push).
Future commits must be checked BEFORE pushing.
