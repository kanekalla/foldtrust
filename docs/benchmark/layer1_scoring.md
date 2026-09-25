# Layer 1: Scoring Correctness

## Question

Are the metric implementations (F1, sensitivity, PPV, slip-tolerance) correct? Do ViennaRNA API calls and parameter loading work as expected? Are the parsers robust?

## Data

Known-answer test cases:
- **Parsing fixtures**: Dot-bracket, bpseq, and CT files under `tests/fixtures/`
- **Pseudoknot example**: Hand-crafted crossing pairs
- **Metric toy cases**: Perfect/zero/partial overlap examples
- **Energy regression**: SARS-CoV-2 FSE (NC_045512.2:13462-13542, 81 nt)
- **Layer 2 sample check**: 600-structure slip-tolerance validation

## Method

Unit tests implemented in `tests/test_benchmark_layer1.py`:

1. **Parsers**:
   - `parse_dotbracket`: `()`, `<>`, `[]`, `{}` nested brackets; unbalanced input raises `StructureParsingError`
   - `parse_bpseq`, `parse_ct`: Load fixtures, verify pair sets match expected values
   - Edge cases: empty structures, single-pair stems, isolated base pairs

2. **Pseudoknot handling**:
   - Greedy nested-only removal: sort pairs by 5' position (tie: longer span), keep if no crossing with already-kept
   - Test case: `{(0,5), (1,6), (2,4), (3,7)}` → keeps `{(0,5), (1,6), (3,7)}`, removes `(2,4)`

3. **Exact vs slip-tolerant matching**:
   - **Exact**: Pair (i,j) correct if (i,j) in reference
   - **Slip-tolerant (±1 on one side)**: Predicted (i,j) matches if any of (i,j), (i±1,j), (i,j±1) in reference; reference pair recovered if any of (i,j), (i±1,j), (i,j±1) predicted
   - Metrics: sensitivity, PPV, F1, MCC
   - **Invariant**: Slip F1 ≥ exact F1 (tested on toy cases and the full Layer 2 sample)

4. **Energy regression** (ViennaRNA 2.7.x parameter loading):
   - Fold FSE with Turner2004 / Andronescu2007 / Langdon2018 using a fresh `RNA.md()` after each `params_load_RNA_*()` call
   - **Expected MFE**: -26.00 / -22.26 / -24.70 kcal/mol
   - Verify `fc.eval_structure(mfe_structure)` equals MFE energy

5. **Base-pair probability checks**:
   - BPP matrix is symmetric: `p[i][j] == p[j][i]`
   - Row sums ≤ 1 + 1e-9 (allowing numerical tolerance)
   - Unpaired probability: `1 - Σ_j p_ij` (summing both triangle halves) ≥ 0

6. **FIRM stem check**:
   - GC hairpin `GCGGGCCC` → all stem pairs have p ≥ 0.85

**Command**:
```bash
pytest tests/test_benchmark_layer1.py -v
```

## Results

All Layer 1 tests pass (57 passed, 3 skipped in full suite as of 2026-09-25).

**Energy regression values** (ViennaRNA 2.7.2, FSE NC_045512.2:13462-13542):
- Turner2004: MFE = -26.00 kcal/mol ✓
- Andronescu2007: MFE = -22.26 kcal/mol ✓
- Langdon2018: MFE = -24.70 kcal/mol ✓

**Slip-tolerance check** (Layer 2 sample, 1800 rows = 600 structures × 3 methods):
- All rows satisfy slip F1 ≥ exact F1 ✓
- Previous bug (slip < exact in 1141/1800 rows) fixed in commit c08de41

**Parser coverage**:
- Dot-bracket: `()`, `<>`, `[]`, `{}`, mixed nesting, unbalanced detection ✓
- bpseq/CT: Fixtures load correctly, pair sets match expected ✓
- Pseudoknot removal: Greedy selection produces expected nested subset ✓

## How to Run

```bash
pytest tests/test_benchmark_layer1.py -v
```

Expected output: 57 passed, 3 skipped (skipped tests are for future ML model checks).

## Limitations

1. **Toy cases only**: Real-structure validation is in Layer 2
2. **ViennaRNA only**: LinearPartition, CONTRAfold, MXfold2 not tested
3. **Turner2004 primary**: Energy regression tests three parameter sets, but Layer 2 accuracy uses Turner2004 only
