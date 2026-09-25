# Layer 1: Scoring Correctness Tests

Unit tests validating metric implementations and ViennaRNA API behavior.

## Command

```bash
pytest tests/test_benchmark_layer1.py -v
```

Or via the benchmark runner:

```bash
python3 scripts/run_layers_1_2_3.py
```

## Outputs

- `layer1_tests.json`: Per-test results with assertions and values
- `layer1_tests.csv`: Same in tabular form

## Tests Included

1. Dot-bracket parsing (nested brackets, unbalanced detection)
2. bpseq/CT format parsing
3. Pseudoknot removal (greedy nested-only selection)
4. Metrics on toy cases (perfect/no/half overlap)
5. Slip ≥ exact assertion
6. Energy regression (FSE sequence, three parameter sets)
7. BPP symmetry and row-sum constraints
8. Unpaired probability validation
9. FIRM tier check (GC hairpin)

See `docs/benchmark/layer1_scoring.md` for details.
