# Vienna.py Subprocess Decision for Layer 5

## Decision: KEEP subprocess approach

## Rationale

The current implementation in `src/foldtrust/vienna.py` uses subprocess for both `fold_with_params()` and `compute_pair_probs_with_params()`. This was implemented to work around ViennaRNA's persistent global parameter state.

**Why subprocess is necessary:**
- ViennaRNA maintains global state for energy parameters that does not fully reset between calls in the same process
- Even after calling `RNA.params_load_RNA_Andronescu2007()` followed by `RNA.md()`, subsequent folds may still use stale parameter state
- The subprocess approach guarantees parameter isolation by running each fold in a fresh Python interpreter

**Evidence it works correctly:**
- `tests/test_param_energies.py` passes all tests
- SARS-CoV-2 FSE (NC_045512.2:13462-13542) correctly produces:
  - Turner2004: -26.00 kcal/mol ✓
  - Andronescu2007: -22.26 kcal/mol ✓
  - Langdon2018: -24.70 kcal/mol ✓
- All three values differ as expected, proving parameter switching works

**Performance consideration:**
Layer 5 requires many folds across conditions (3 parameter sets × 4 temperatures × 5 cases, plus window context sweeps). The subprocess overhead (~100-200ms per fold) is acceptable for a benchmark that runs in minutes.

**Alternative considered:**
Revert to in-process folding with careful `RNA.md()` creation after each `params_load`. However, this proved unreliable in previous implementations (see git history), and the working subprocess approach should not be changed without a compelling reason.

## Conclusion

Keep the subprocess approach as-is. It is reliable, tested, and produces correct results.
