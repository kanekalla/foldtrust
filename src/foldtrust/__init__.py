"""
FoldTrust: Ensemble-first RNA structure reliability reports

A scanpy-style toolkit for RNA secondary structure prediction with ensemble uncertainty quantification.

Import convention::

    import foldtrust as ft

Core data structure::

    fd = ft.io.read_fasta("sequence.fa")  # FoldData object

Workflow::

    # I/O
    fd = ft.io.read_fasta("sequence.fa")
    fd = ft.io.read_case_directory("data/cases/sars2-fse")
    collection = ft.io.read_case_batch("data/cases")
    
    # Analysis tools
    ft.tl.fold_mfe(fd)
    ft.tl.compute_ensemble(fd)
    ft.tl.call_tiers(fd)
    ft.tl.compute_unpaired_probs(fd)
    # Or run complete pipeline:
    fd = ft.tl.run_pipeline(fd, temperature=37.0)
    
    # Plotting
    ft.pl.arc_plot(fd, save="arc.png")
    ft.pl.heatmap(fd, save="heatmap.png")
    ft.pl.tier_distribution(fd)
    ft.pl.shape_track(fd)  # if SHAPE data loaded
    ft.pl.disease_summary(collection)

Data persistence::

    fd.write_h5("output.h5")
    fd2 = ft.FoldData.read_h5("output.h5")

"""

__version__ = "0.2.0"

# Core data structures
from foldtrust._core import FoldData, FoldDataCollection

# Namespace imports
from foldtrust import io
from foldtrust import tl
from foldtrust import pl

# Legacy CLI support
from foldtrust.core import process_sequence
from foldtrust.utils import read_fasta as _legacy_read_fasta, find_case_directories

__all__ = [
    # Core
    "FoldData",
    "FoldDataCollection",
    # Namespaces
    "io",
    "tl",
    "pl",
    # Version
    "__version__",
]
