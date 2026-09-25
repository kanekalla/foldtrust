# FoldTrust API Documentation

**Scanpy-style API for RNA structure ensemble analysis**

## Installation

```bash
pip install -e .
pip install h5py pyyaml zarr
```

## Quick Start

```python
import foldtrust as ft

# Load a sequence
fd = ft.io.read_fasta("data/cases/sars2-fse/sequence.fa")

# Run complete analysis pipeline
fd = ft.tl.run_pipeline(fd, temperature=37.0)

# Access results
print(fd.structures['mfe'])  # MFE structure
print(fd.obs['tier'])  # Per-nucleotide tier classifications
print(fd.obs['unpaired_prob'])  # Unpaired probabilities

# Plot results
ft.pl.arc_plot(fd, save="arc.png")
ft.pl.heatmap(fd, save="heatmap.png")
ft.pl.tier_distribution(fd)

# Save to disk
fd.write_h5("output.h5")
fd2 = ft.FoldData.read_h5("output.h5")
```

## API Structure

### Core Data Container

**`FoldData`** - AnnData-style container for RNA structure analysis

- `sequence`: RNA sequence (str)
- `name`: Identifier (str)
- `obs`: Per-nucleotide annotations (pd.DataFrame)
  - `tier`: FIRM/SOFT/FLOPPY classification
  - `unpaired_prob`: Unpaired probability
  - `pair_prob`: Pairing probability (if paired)
  - `shape_*`: SHAPE reactivity data
- `obsp`: Pairwise data (dict of np.ndarray)
  - `pair_probs`: Pair probability matrix (n × n)
- `structures`: Named secondary structures (dict)
  - `mfe`: Minimum free energy structure
  - `reference`: Known reference structure (if available)
- `layers`: Named data layers for different conditions (dict)
- `uns`: Unstructured metadata (dict)

**`FoldDataCollection`** - Collection of FoldData objects for batch analysis

### Namespaces

#### `ft.io` - Input/Output

```python
ft.io.read_fasta(filename)  # Load FASTA file
ft.io.read_case_directory(case_dir)  # Load case with metadata
ft.io.read_case_batch(cases_dir)  # Load multiple cases
ft.io.read_reference_structure(seq_file, struct_file)  # Load with reference
ft.io.load_shape_data(fd, shape_file, start, end)  # Add SHAPE data
```

#### `ft.tl` - Analysis Tools

```python
ft.tl.fold_mfe(fd)  # Predict MFE structure
ft.tl.compute_ensemble(fd)  # Compute pair probabilities
ft.tl.call_tiers(fd)  # Classify tiers (FIRM/SOFT/FLOPPY)
ft.tl.compute_unpaired_probs(fd)  # Compute unpaired probabilities
ft.tl.run_pipeline(fd)  # Run complete pipeline
ft.tl.run_pipeline_batch(collection)  # Batch processing
```

#### `ft.pl` - Plotting

```python
ft.pl.arc_plot(fd)  # Arc diagram of structure
ft.pl.heatmap(fd)  # Pair probability heatmap
ft.pl.tier_distribution(fd)  # Tier composition bar chart
ft.pl.shape_track(fd)  # SHAPE reactivity vs unpaired probability
ft.pl.disease_summary(collection)  # Summary across multiple windows
```

## Data Persistence

### HDF5 (recommended)

```python
fd.write_h5("output.h5")
fd = ft.FoldData.read_h5("output.h5")
```

### Zarr (alternative)

```python
fd.write_zarr("output.zarr")
fd = ft.FoldData.read_zarr("output.zarr")
```

### Collections

```python
collection.write_h5("output_dir/")
collection = ft.FoldDataCollection.read_h5("output_dir/")
```

## Examples

### Disease Window Analysis

```python
import foldtrust as ft

# Load all disease cases
collection = ft.io.read_case_batch("data/cases")

# Analyze each case
collection = ft.tl.run_pipeline_batch(collection, temperature=37.0)

# Visualize summary
ft.pl.disease_summary(collection, save="disease_summary.png")

# Access individual results
for name, fd in collection:
    print(f"{name}: {fd.uns['n_firm']} FIRM stems")
```

### SHAPE Data Integration

```python
import foldtrust as ft

# Load SARS-CoV-2 FSE
fd = ft.io.read_case_directory("data/cases/sars2-fse")

# Run analysis
fd = ft.tl.run_pipeline(fd)

# Load SHAPE data
fd = ft.io.load_shape_data(
    fd,
    "/tmp/SARS_CoV-2_shape_comparison/SHAPE data/zhang_invivo_reactivity.csv",
    start=13468,
    end=13638,
    dataset_name="shape_zhang"
)

# Compare SHAPE with unpaired probabilities
ft.pl.shape_track(fd, shape_key="shape_zhang", save="shape_track.png")

# Compute correlation
from scipy.stats import spearmanr
valid = ~fd.obs['shape_zhang'].isna()
rho, pval = spearmanr(
    fd.obs.loc[valid, 'unpaired_prob'],
    fd.obs.loc[valid, 'shape_zhang']
)
print(f"Spearman ρ = {rho:.3f}, p = {pval:.2e}")
```

### Reference Structure Benchmark

```python
import foldtrust as ft

# Load reference structure
fd = ft.io.read_reference_structure(
    "benchmarks/reference_data/5S_rRNA_ecoli.fa",
    "benchmarks/reference_data/5S_rRNA_ecoli.ct"
)

# Predict structure
fd = ft.tl.run_pipeline(fd)

# Compare MFE to reference
def compare_structures(ref_struct, pred_struct):
    """Compare two structures and return metrics."""
    # Parse pairs from structures
    ref_pairs = set()
    stack = []
    for i, c in enumerate(ref_struct):
        if c == '(':
            stack.append(i)
        elif c == ')' and stack:
            j = stack.pop()
            ref_pairs.add((j, i))
    
    pred_pairs = set()
    stack = []
    for i, c in enumerate(pred_struct):
        if c == '(':
            stack.append(i)
        elif c == ')' and stack:
            j = stack.pop()
            pred_pairs.add((j, i))
    
    tp = len(ref_pairs & pred_pairs)
    fp = len(pred_pairs - ref_pairs)
    fn = len(ref_pairs - pred_pairs)
    
    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0
    ppv = tp / (tp + fp) if (tp + fp) > 0 else 0
    f1 = 2 * tp / (2 * tp + fp + fn) if (2 * tp + fp + fn) > 0 else 0
    
    return {"sensitivity": sensitivity, "ppv": ppv, "f1": f1}

metrics = compare_structures(fd.structures['reference'], fd.structures['mfe'])
print(f"F1 score: {metrics['f1']:.3f}")
```

## CLI Compatibility

The CLI commands (`report`, `batch`, `demo`, `benchmark`) remain fully functional and now use the new API internally.

```bash
foldtrust report data/cases/sars2-fse/sequence.fa -o output/
foldtrust batch data/cases/ -o output/batch/
foldtrust demo
foldtrust benchmark all -o benchmarks/outputs/
```

## Development

### Running Tests

```python
pytest tests/
```

### Adding New Analysis Tools

Follow the scanpy pattern:

```python
# In src/foldtrust/tl.py

def my_new_tool(fd: FoldData, param: float = 1.0, in_place: bool = True) -> Optional[FoldData]:
    """
    My new analysis tool.
    
    Parameters
    ----------
    fd : FoldData
        Input data
    param : float
        Analysis parameter
    in_place : bool
        Modify in place or return copy
    """
    if not in_place:
        fd = fd.copy()
    
    # Do analysis, write results to fd.obs, fd.obsp, or fd.uns
    fd.obs['my_result'] = ...
    
    if not in_place:
        return fd
```

## References

- **AnnData/scanpy design pattern**: Wolf et al., Genome Biol 2018. doi:10.1186/s13059-017-1382-0
- **ViennaRNA Package**: Lorenz et al., Algorithms Mol Biol 2011. doi:10.1186/1748-7188-6-26
