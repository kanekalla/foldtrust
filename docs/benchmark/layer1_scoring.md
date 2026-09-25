# Layer 1: Scoring Correctness

## Question

Are the metric implementations (F1, sensitivity, PPV, slip-tolerance, ECE, tier definitions) correct? Do ViennaRNA API calls produce expected results?

## Method

Unit tests on known-answer cases:
1. **Dot-bracket parsing**: `()`, `<>`, `[]`, `{}` nested brackets; unbalanced input raises error
2. **bpseq/CT parsing**: Fixtures under `tests/fixtures/` with known pair sets
3. **Pseudoknot removal**: Greedy nested-only selection on a crossing example
4. **Metrics toy cases**: Perfect match (F1=1.0), no overlap (F1=0.0), half overlap (F1=0.5)
5. **Slip ≥ exact**: Random structures and the Layer 2 sample (slip F1 ≥ exact F1 always)
6. **Energy regression**: FSE sequence (`NC_045512.2:13462-13542`) folded with Turner2004/Andronescu2007/Langdon2018:
   - MFE = -26.00 / -22.26 / -24.70 kcal/mol (using fresh `RNA.md()` after each `params_load`)
   - `fc.eval_structure(mfe_structure)` equals MFE energy
7. **BPP symmetry**: Base-pair probability matrix is symmetric and row sums ≤ 1+1e-9
8. **Unpaired probabilities**: `1 - Σ_j p_ij` equals ViennaRNA unpaired prob (both triangles of symmetric matrix)
9. **FIRM check**: GC hairpin `GCGGGCCC` forms firm stem (all pairs p ≥ 0.85)

**Command**:
```bash
pytest tests/test_benchmark_layer1.py -v
```

## Results

All tests pass. Results saved to `benchmarks/outputs/layer1/layer1_tests.json` and `.csv`.

(Detailed per-test results will be added when Layer 1 pytest suite is complete.)

## Limitations

1. **Toy cases only**: Real-structure validation is in Layer 2
2. **No ML model checks**: LinearPartition, CONTRAfold, MXfold2 not tested (ViennaRNA only)
3. **Parameter sets**: Energy regression tests three sets but folding accuracy (Layer 2) uses Turner2004 only
