"""ViennaRNA integration for RNA folding."""

import json
import re
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

try:
    import RNA
    HAS_RNA_PYTHON = True
except ImportError:
    HAS_RNA_PYTHON = False


class ViennaRNAError(Exception):
    """Raised when ViennaRNA calls fail."""

    pass


def check_viennarna() -> bool:
    """Check if ViennaRNA is available."""
    try:
        result = subprocess.run(["RNAfold", "--version"], capture_output=True, text=True, timeout=5)
        return result.returncode == 0
    except (subprocess.SubprocessError, FileNotFoundError):
        return False


def fold_mfe(sequence: str) -> Tuple[str, float]:
    """
    Compute MFE structure using RNAfold.

    Returns:
        (structure, energy) tuple where structure is dot-bracket notation
    """
    if not check_viennarna():
        raise ViennaRNAError("RNAfold not found. Install ViennaRNA: brew install viennarna")

    try:
        result = subprocess.run(
            ["RNAfold", "--noPS"],
            input=f">seq\n{sequence}\n",
            capture_output=True,
            text=True,
            timeout=60,
        )

        if result.returncode != 0:
            raise ViennaRNAError(f"RNAfold failed: {result.stderr}")

        lines = result.stdout.strip().split("\n")
        structure_line = lines[-1]

        match = re.search(r"([.()\[\]{}]+)\s+\(\s*(-?\d+\.\d+)\s*\)", structure_line)
        if not match:
            raise ViennaRNAError(f"Could not parse RNAfold output: {structure_line}")

        structure = match.group(1)
        energy = float(match.group(2))

        return structure, energy

    except subprocess.TimeoutExpired:
        raise ViennaRNAError("RNAfold timed out")
    except Exception as e:
        raise ViennaRNAError(f"RNAfold error: {e}")


def compute_pair_probabilities(sequence: str) -> np.ndarray:
    """
    Compute base-pair probability matrix using RNAfold partition function.

    Returns:
        NxN matrix where element [i,j] is probability of pairing between positions i and j
    """
    try:
        import RNA
    except ImportError:
        raise ViennaRNAError("ViennaRNA Python package not found. Install: pip install ViennaRNA")

    n = len(sequence)
    prob_matrix = np.zeros((n, n))

    try:
        # Create fold compound
        md = RNA.md()
        md.uniq_ML = 1
        fc = RNA.fold_compound(sequence, md)
        
        # Compute partition function
        fc.pf()
        
        # Get base-pair probabilities
        bpp = fc.bpp()
        
        # Convert to matrix (ViennaRNA uses 1-based indexing)
        # bpp[i][j] is upper-triangular (only i < j has values)
        # Make it symmetric by copying to both (i,j) and (j,i)
        for i in range(1, n + 1):
            for j in range(i + 1, n + 1):
                prob = bpp[i][j]
                if prob > 0:
                    prob_matrix[i - 1, j - 1] = prob
                    prob_matrix[j - 1, i - 1] = prob

        return prob_matrix

    except Exception as e:
        raise ViennaRNAError(f"Partition function error: {e}")


def parse_stems(structure: str, prob_matrix: np.ndarray, min_stem_length: int = 2) -> List[Dict]:
    """
    Parse stems from MFE structure and compute reliability metrics.

    Args:
        structure: Dot-bracket structure
        prob_matrix: Base-pair probability matrix
        min_stem_length: Minimum number of consecutive base pairs to count as stem

    Returns:
        List of stem dictionaries with positions, length, mean probability, and flag
    """
    len(structure)
    pairs = []
    stack = []

    for i, char in enumerate(structure):
        if char in "([{":
            stack.append(i)
        elif char in ")]}":
            if stack:
                j = stack.pop()
                pairs.append((j, i))

    pairs.sort()

    stems = []
    current_stem = []

    for k, (i, j) in enumerate(pairs):
        if not current_stem:
            current_stem.append((i, j))
        else:
            last_i, last_j = current_stem[-1]
            if i == last_i + 1 and j == last_j - 1:
                current_stem.append((i, j))
            else:
                if len(current_stem) >= min_stem_length:
                    stems.append(current_stem)
                current_stem = [(i, j)]

    if len(current_stem) >= min_stem_length:
        stems.append(current_stem)

    stem_data = []
    for stem_id, stem_pairs in enumerate(stems, 1):
        probs = [prob_matrix[i, j] for i, j in stem_pairs]
        mean_prob = np.mean(probs)

        if mean_prob >= 0.85:
            flag = "firm"
        elif mean_prob >= 0.5:
            flag = "soft"
        else:
            flag = "floppy"

        stem_data.append(
            {
                "id": stem_id,
                "pairs": stem_pairs,
                "length": len(stem_pairs),
                "mean_prob": mean_prob,
                "flag": flag,
                "positions": f"{stem_pairs[0][0]+1}-{stem_pairs[0][1]+1} ... {stem_pairs[-1][0]+1}-{stem_pairs[-1][1]+1}",
            }
        )

    return stem_data


def fold_with_params(
    sequence: str,
    param_set: str = "Turner2004",
    temperature: float = 37.0,
) -> Tuple[str, float]:
    """
    Fold RNA with specific parameter set and temperature using Python API.
    
    ViennaRNA has persistent global state for parameters that doesn't fully
    reset between calls in the same process. To get reliable parameter switching,
    we run each fold in a fresh subprocess.
    
    Args:
        sequence: RNA sequence
        param_set: One of "Turner2004", "Andronescu2007", "Langdon2018"
        temperature: Temperature in Celsius
    
    Returns:
        (structure, mfe_energy) tuple
        
    Raises:
        ViennaRNAError: If Python API not available or parameters invalid
    """
    if not HAS_RNA_PYTHON:
        raise ViennaRNAError(
            "ViennaRNA Python bindings not available. "
            "Install with: pip install ViennaRNA"
        )
    
    if param_set not in ["Turner2004", "Andronescu2007", "Langdon2018"]:
        raise ViennaRNAError(f"Unknown parameter set: {param_set}")
    
    # Run in subprocess to avoid parameter state issues
    import sys
    code = f"""
import RNA
import sys

sequence = {sequence!r}
param_set = {param_set!r}
temperature = {temperature}

# Load parameters
if param_set == "Turner2004":
    RNA.params_load_RNA_Turner2004()
elif param_set == "Andronescu2007":
    RNA.params_load_RNA_Andronescu2007()
elif param_set == "Langdon2018":
    RNA.params_load_RNA_Langdon2018()

# Create fold compound with fresh md
md = RNA.md()
md.temperature = temperature
fc = RNA.fold_compound(sequence, md)

# Compute MFE
structure, mfe = fc.mfe()

print(f"{{structure}}|{{mfe}}")
"""
    
    try:
        result = subprocess.run(
            [sys.executable, "-c", code],
            capture_output=True,
            text=True,
            timeout=30,
            check=True,
        )
        
        output = result.stdout.strip()
        structure, mfe_str = output.split("|")
        mfe = float(mfe_str)
        
        return structure, mfe
        
    except subprocess.CalledProcessError as e:
        raise ViennaRNAError(f"Fold failed: {e.stderr}")
    except subprocess.TimeoutExpired:
        raise ViennaRNAError("Fold timed out")
    except Exception as e:
        raise ViennaRNAError(f"Fold error: {e}")


def compute_pair_probs_with_params(
    sequence: str,
    param_set: str = "Turner2004", 
    temperature: float = 37.0,
) -> np.ndarray:
    """
    Compute base-pair probabilities with specific parameters.
    
    Uses subprocess to avoid ViennaRNA parameter state issues.
    
    Args:
        sequence: RNA sequence
        param_set: One of "Turner2004", "Andronescu2007", "Langdon2018"
        temperature: Temperature in Celsius
        
    Returns:
        NxN probability matrix
        
    Raises:
        ViennaRNAError: If Python API not available
    """
    if not HAS_RNA_PYTHON:
        raise ViennaRNAError(
            "ViennaRNA Python bindings not available. "
            "Install with: pip install ViennaRNA"
        )
    
    if param_set not in ["Turner2004", "Andronescu2007", "Langdon2018"]:
        raise ViennaRNAError(f"Unknown parameter set: {param_set}")
    
    # Run in subprocess
    import sys
    import json
    code = f"""
import RNA
import json

sequence = {sequence!r}
param_set = {param_set!r}
temperature = {temperature}

# Load parameters
if param_set == "Turner2004":
    RNA.params_load_RNA_Turner2004()
elif param_set == "Andronescu2007":
    RNA.params_load_RNA_Andronescu2007()
elif param_set == "Langdon2018":
    RNA.params_load_RNA_Langdon2018()

# Create fold compound with fresh md
md = RNA.md()
md.temperature = temperature
fc = RNA.fold_compound(sequence, md)

# Compute partition function
fc.pf()

# Extract base-pair probabilities
n = len(sequence)
bpp = fc.bpp()

# Build sparse pairs list (to avoid huge JSON)
pairs = []
for i in range(1, n + 1):
    for j in range(i + 1, n + 1):
        prob = bpp[i][j]
        if prob > 0:
            pairs.append([i-1, j-1, prob])

print(json.dumps({{"n": n, "pairs": pairs}}))
"""
    
    try:
        result = subprocess.run(
            [sys.executable, "-c", code],
            capture_output=True,
            text=True,
            timeout=60,
            check=True,
        )
        
        data = json.loads(result.stdout)
        n = data["n"]
        pairs = data["pairs"]
        
        # Reconstruct matrix
        prob_matrix = np.zeros((n, n))
        for i, j, prob in pairs:
            prob_matrix[i, j] = prob
            prob_matrix[j, i] = prob
        
        return prob_matrix
        
    except subprocess.CalledProcessError as e:
        raise ViennaRNAError(f"Partition function failed: {e.stderr}")
    except subprocess.TimeoutExpired:
        raise ViennaRNAError("Partition function timed out")
    except Exception as e:
        raise ViennaRNAError(f"Partition function error: {e}")


def load_param_file(param_file: Path) -> None:
    """
    Load ViennaRNA parameter file.
    
    Args:
        param_file: Path to .par file
        
    Raises:
        ViennaRNAError: If Python API not available or file not found
    """
    if not HAS_RNA_PYTHON:
        raise ViennaRNAError(
            "ViennaRNA Python bindings not available. "
            "Install with: pip install ViennaRNA"
        )
    
    if not param_file.exists():
        raise ViennaRNAError(f"Parameter file not found: {param_file}")
    
    RNA.params_load(str(param_file))
