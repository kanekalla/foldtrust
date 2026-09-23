"""ViennaRNA integration for RNA folding."""

import re
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np


class ViennaRNAError(Exception):
    """Raised when ViennaRNA calls fail."""
    pass


def check_viennarna() -> bool:
    """Check if ViennaRNA is available."""
    try:
        result = subprocess.run(
            ["RNAfold", "--version"],
            capture_output=True,
            text=True,
            timeout=5
        )
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
            timeout=60
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
    if not check_viennarna():
        raise ViennaRNAError("RNAfold not found. Install ViennaRNA: brew install viennarna")
    
    n = len(sequence)
    prob_matrix = np.zeros((n, n))
    
    try:
        result = subprocess.run(
            ["RNAfold", "-p", "--noPS"],
            input=f">seq\n{sequence}\n",
            capture_output=True,
            text=True,
            timeout=120
        )
        
        if result.returncode != 0:
            raise ViennaRNAError(f"RNAfold partition function failed: {result.stderr}")
        
        dp_file = Path("seq_dp.ps")
        if not dp_file.exists():
            raise ViennaRNAError("RNAfold did not generate dot plot PostScript file")
        
        with open(dp_file, 'r') as f:
            in_data_section = False
            for line in f:
                line = line.strip()
                
                if line == "/sequence { (":
                    in_data_section = False
                elif line.startswith("%start of base pair probability data"):
                    in_data_section = True
                    continue
                elif in_data_section:
                    if line.startswith("showpage"):
                        break
                    
                    parts = line.split()
                    if len(parts) >= 4 and parts[3] == "ubox":
                        try:
                            i = int(parts[0]) - 1
                            j = int(parts[1]) - 1
                            sqrt_prob = float(parts[2])
                            prob = sqrt_prob * sqrt_prob
                            
                            if 0 <= i < n and 0 <= j < n:
                                prob_matrix[i, j] = prob
                                prob_matrix[j, i] = prob
                        except (ValueError, IndexError):
                            continue
        
        dp_file.unlink()
        
        dot_file = Path("dot.ps")
        if dot_file.exists():
            dot_file.unlink()
        
        return prob_matrix
        
    except subprocess.TimeoutExpired:
        raise ViennaRNAError("RNAfold partition function timed out")
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
    n = len(structure)
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
        
        stem_data.append({
            "id": stem_id,
            "pairs": stem_pairs,
            "length": len(stem_pairs),
            "mean_prob": mean_prob,
            "flag": flag,
            "positions": f"{stem_pairs[0][0]+1}-{stem_pairs[0][1]+1} ... {stem_pairs[-1][0]+1}-{stem_pairs[-1][1]+1}"
        })
    
    return stem_data
