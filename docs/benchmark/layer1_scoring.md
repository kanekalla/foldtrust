# Layer 1: Scoring Correctness

## Question

Do our base-pair scoring functions work correctly? Can they handle bracket notations, bpseq, and CT formats? Do tier assignments behave as expected for strong vs weak structures?

## Method

### Parser Tests

We test three structure format parsers:

1. **Dot-bracket notation**: Parse `()`, `<>`, `[]`, `{}` as nested pairs following WUSS conventions
2. **bpseq format**: `index nucleotide pair_index` (1-indexed, 0 = unpaired)
3. **CT format**: Connectivity table with pairing columns

Test cases include:
- Simple hairpins: `(((...)))..(((...)))`
- Nested brackets: `(((<...>)))`
- Multiple bracket types: `(((...))).[[...]].{{...}}`

### Metrics Tests

We verify correctness of sensitivity, PPV, F1, and MCC computation with known test cases:

- **Perfect prediction**: All metrics = 1.0
- **Mixed errors**: Known TP, FP, FN counts
- **Slip tolerance**: +/-1 nucleotide slippage increases TP count

### Tier Regression Tests

We test that a **strong GC hairpin** (`GCGCGCAAAAGCGCGC`) produces FIRM pairs (p ≥ 0.85):

- Compute ViennaRNA MFE and partition function
- Extract MFE pairs and their probabilities
- Verify that most pairs are classified as FIRM

**Rationale**: A 6-bp GC stem should be thermodynamically stable with high pair probabilities in the ensemble.

**Note on removed references**: The old benchmark contained a hand-edited tRNA-Phe reference file that was altered to match a prediction. This file has been **removed** from the benchmark. All reference structures in Layers 2 and 3 come verbatim from the cited sources (ArchiveII, bpRNA, Rfam, SPOT-RNA) and are never edited toward predictions.

## Data

No external data required for Layer 1 tests (synthetic test cases only).

## Command

```bash
cd /workspace
python3 src/foldtrust/benchmark/layer1_scoring.py
```

Or via the main runner:

```bash
python3 scripts/run_layers_1_2_3.py
```

## Results

All tests passed:

- ✓ Dot-bracket parser
- ✓ Nested bracket parser
- ✓ Multiple bracket types
- ✓ Perfect prediction test
- ✓ FP/FN test
- ✓ Slip tolerance test
- ✓ Tier regression test (strong GC hairpin: 6 pairs, 6 FIRM)

**Strong GC hairpin results**:
- Sequence: `GCGCGCAAAAGCGCGC`
- MFE structure: `((((((....)))))))`
- MFE energy: -10.90 kcal/mol
- All 6 stem pairs have p ≥ 0.85 (FIRM)

Detailed results: `benchmarks/outputs/layer1/layer1_tests.json`

## Interpretation

The scoring functions behave correctly:

1. **Parsers** handle all tested notation formats correctly
2. **Metrics** compute TP/FP/FN/TN correctly and derive sensitivity/PPV/F1/MCC
3. **Slip tolerance** correctly allows +/-1 slippage
4. **Tier regression** confirms that thermodynamically strong structures produce FIRM pairs

These tests provide confidence that Layers 2 and 3 results are not artifacts of scoring bugs.

## Limitations

- **Tier test uses a single strong hairpin**: More diverse test structures (internal loops, bulges, multi-branch loops) would strengthen the regression suite
- **No weak structure test**: We did not test that a known floppy structure produces FLOPPY pairs (could add in future)
- **Pseudoknot handling not tested**: Layers 2 and 3 remove pseudoknots; Layer 1 does not explicitly test this removal

## Files Generated

```
benchmarks/outputs/layer1/
└── layer1_tests.json      # Test results summary
```
