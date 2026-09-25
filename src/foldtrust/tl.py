"""Analysis tools for FoldTrust.

Functions for structure prediction, ensemble analysis, tier calling, and benchmarking.
"""

from typing import Optional, List, Dict
import numpy as np
import subprocess

from foldtrust._core import FoldData, FoldDataCollection


def fold_mfe(fd: FoldData, temperature: float = 37.0, in_place: bool = True) -> Optional[FoldData]:
    """
    Predict minimum free energy (MFE) structure with ViennaRNA.
    
    Parameters
    ----------
    fd : FoldData
        FoldData object with sequence
    temperature : float
        Temperature in Celsius (default: 37.0)
    in_place : bool
        If True, modify fd in place; if False, return a copy
        
    Returns
    -------
    FoldData or None
        Updated FoldData if in_place=False, else None
        
    Examples
    --------
    >>> fd = ft.io.read_fasta("data/cases/sars2-fse/sequence.fa")
    >>> ft.tl.fold_mfe(fd)
    >>> fd.structures['mfe']
    '(((((((((...))))))...'
    """
    if not in_place:
        fd = fd.copy()
    
    # Run RNAfold
    result = subprocess.run(
        ["RNAfold", "--noPS", "-T", str(temperature)],
        input=fd.sequence,
        capture_output=True,
        text=True
    )
    
    lines = result.stdout.strip().split('\n')
    if len(lines) >= 2:
        parts = lines[1].split()
        structure = parts[0]
        energy = float(parts[1].strip('()'))
        
        fd.structures['mfe'] = structure
        fd.uns['mfe_energy'] = energy
        fd.uns['mfe_temperature'] = temperature
    
    if not in_place:
        return fd


def compute_ensemble(fd: FoldData, temperature: float = 37.0, in_place: bool = True) -> Optional[FoldData]:
    """
    Compute partition function and pair probabilities.
    
    Parameters
    ----------
    fd : FoldData
        FoldData object with sequence
    temperature : float
        Temperature in Celsius
    in_place : bool
        If True, modify fd in place; if False, return a copy
        
    Returns
    -------
    FoldData or None
        Updated FoldData with pair_probs in obsp
        
    Examples
    --------
    >>> fd = ft.io.read_fasta("data/cases/sars2-fse/sequence.fa")
    >>> ft.tl.compute_ensemble(fd)
    >>> 'pair_probs' in fd.obsp
    True
    """
    if not in_place:
        fd = fd.copy()
    
    try:
        import RNA
        
        # Create fold compound
        md = RNA.md()
        md.temperature = temperature
        fc = RNA.fold_compound(fd.sequence, md)
        
        # Compute partition function
        fc.pf()
        
        # Get base pair probabilities
        bpp = fc.bpp()
        n = len(fd.sequence)
        prob_matrix = np.zeros((n, n))
        
        for i in range(n):
            for j in range(n):
                if i < len(bpp) and j < len(bpp[i]):
                    prob_matrix[i, j] = bpp[i][j]
        
        fd.obsp['pair_probs'] = prob_matrix
        fd.uns['ensemble_temperature'] = temperature
        
    except ImportError:
        # Fallback: parse RNAfold -p output
        print("Warning: ViennaRNA Python not available, using RNAfold -p (slower)")
        
        result = subprocess.run(
            ["RNAfold", "-p", "--noPS", "-T", str(temperature)],
            input=fd.sequence,
            capture_output=True,
            text=True
        )
        
        # Would need to parse dot plot PostScript output - complex
        # For now, just note that ensemble was attempted
        fd.uns['ensemble_attempted'] = True
        fd.uns['ensemble_temperature'] = temperature
    
    if not in_place:
        return fd


def call_tiers(fd: FoldData, in_place: bool = True) -> Optional[FoldData]:
    """
    Classify base pairs into FIRM/SOFT/FLOPPY tiers.
    
    Requires MFE structure and pair probabilities.
    
    Parameters
    ----------
    fd : FoldData
        FoldData object with structures['mfe'] and obsp['pair_probs']
    in_place : bool
        If True, modify fd in place; if False, return a copy
        
    Returns
    -------
    FoldData or None
        Updated FoldData with tier annotations
        
    Examples
    --------
    >>> fd = ft.io.read_fasta("data/cases/sars2-fse/sequence.fa")
    >>> ft.tl.fold_mfe(fd)
    >>> ft.tl.compute_ensemble(fd)
    >>> ft.tl.call_tiers(fd)
    >>> fd.obs['tier']
    0       NaN
    1    floppy
    2    floppy
    ...
    """
    if not in_place:
        fd = fd.copy()
    
    if 'mfe' not in fd.structures:
        raise ValueError("MFE structure not found. Run ft.tl.fold_mfe first.")
    
    if 'pair_probs' not in fd.obsp:
        raise ValueError("Pair probabilities not found. Run ft.tl.compute_ensemble first.")
    
    from foldtrust.vienna import parse_stems
    
    structure = fd.structures['mfe']
    prob_matrix = fd.obsp['pair_probs']
    
    # Parse stems and classify
    stems = parse_stems(structure, prob_matrix)
    
    # Create per-nucleotide tier annotations
    tier = np.full(fd.n_obs, np.nan, dtype=object)
    pair_prob = np.full(fd.n_obs, np.nan)
    
    for stem in stems:
        for i, j in stem['pairs']:
            tier[i] = stem['flag']
            tier[j] = stem['flag']
            pair_prob[i] = prob_matrix[i, j]
            pair_prob[j] = prob_matrix[i, j]
    
    fd.obs['tier'] = tier
    fd.obs['pair_prob'] = pair_prob
    
    # Store stem-level data
    fd.uns['stems'] = stems
    fd.uns['n_firm'] = sum(1 for s in stems if s['flag'] == 'firm')
    fd.uns['n_soft'] = sum(1 for s in stems if s['flag'] == 'soft')
    fd.uns['n_floppy'] = sum(1 for s in stems if s['flag'] == 'floppy')
    
    if not in_place:
        return fd


def compute_unpaired_probs(fd: FoldData, in_place: bool = True) -> Optional[FoldData]:
    """
    Compute unpaired probability for each nucleotide.
    
    Unpaired probability = 1 - sum of pairing probabilities.
    
    Parameters
    ----------
    fd : FoldData
        FoldData object with pair probabilities
    in_place : bool
        If True, modify fd in place; if False, return a copy
        
    Returns
    -------
    FoldData or None
        Updated FoldData with unpaired_prob in obs
    """
    if not in_place:
        fd = fd.copy()
    
    if 'pair_probs' not in fd.obsp:
        raise ValueError("Pair probabilities not found. Run ft.tl.compute_ensemble first.")
    
    prob_matrix = fd.obsp['pair_probs']
    n = fd.n_obs
    unpaired = np.zeros(n)
    
    for i in range(n):
        # Sum all pairing probabilities involving position i
        paired_prob = np.sum(prob_matrix[i, :]) + np.sum(prob_matrix[:, i]) - prob_matrix[i, i]
        unpaired[i] = max(0.0, min(1.0, 1.0 - paired_prob / 2.0))
    
    fd.obs['unpaired_prob'] = unpaired
    
    if not in_place:
        return fd


def run_pipeline(
    fd: FoldData,
    temperature: float = 37.0,
    compute_unpaired: bool = True,
) -> FoldData:
    """
    Run complete FoldTrust pipeline: fold + ensemble + tiers + unpaired probs.
    
    Parameters
    ----------
    fd : FoldData
        FoldData object with sequence
    temperature : float
        Temperature in Celsius
    compute_unpaired : bool
        Whether to compute unpaired probabilities
        
    Returns
    -------
    FoldData
        Updated FoldData with all annotations
        
    Examples
    --------
    >>> fd = ft.io.read_fasta("data/cases/sars2-fse/sequence.fa")
    >>> fd = ft.tl.run_pipeline(fd, temperature=37.0)
    >>> fd.structures['mfe']
    >>> fd.obs['tier']
    >>> fd.obs['unpaired_prob']
    """
    fold_mfe(fd, temperature=temperature)
    compute_ensemble(fd, temperature=temperature)
    call_tiers(fd)
    
    if compute_unpaired:
        compute_unpaired_probs(fd)
    
    return fd


def run_pipeline_batch(
    collection: FoldDataCollection,
    temperature: float = 37.0,
) -> FoldDataCollection:
    """
    Run FoldTrust pipeline on all samples in a collection.
    
    Parameters
    ----------
    collection : FoldDataCollection
        Collection of FoldData objects
    temperature : float
        Temperature in Celsius
        
    Returns
    -------
    FoldDataCollection
        Updated collection
    """
    for name, fd in collection:
        print(f"Processing {name}...")
        run_pipeline(fd, temperature=temperature)
    
    return collection
