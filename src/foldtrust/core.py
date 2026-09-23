"""Core processing logic for FoldTrust."""

from pathlib import Path
from typing import Dict

import yaml

from foldtrust.utils import read_fasta, compute_gc_content
from foldtrust.vienna import fold_mfe, compute_pair_probabilities, parse_stems, check_viennarna
from foldtrust.report import generate_html_report, generate_markdown_report
from foldtrust.visualization import create_heatmap


def process_sequence(
    sequence_path: Path,
    output_dir: Path,
    generate_html: bool = True,
    generate_markdown: bool = True
) -> Dict:
    """
    Process a sequence file and generate structure reliability report.
    
    Args:
        sequence_path: Path to FASTA file
        output_dir: Directory for output files
        generate_html: Whether to generate HTML report
        generate_markdown: Whether to generate Markdown report
    
    Returns:
        Dictionary with analysis results
    """
    if not check_viennarna():
        raise RuntimeError(
            "ViennaRNA not found. Install it:\n"
            "  macOS:  brew install viennarna\n"
            "  Ubuntu: sudo apt-get install vienna-rna"
        )
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    header, sequence = read_fasta(sequence_path)
    
    meta_path = sequence_path.parent / "meta.yaml"
    metadata = {}
    if meta_path.exists():
        with open(meta_path, 'r') as f:
            metadata = yaml.safe_load(f) or {}
    
    gc_content = compute_gc_content(sequence)
    
    structure, mfe_energy = fold_mfe(sequence)
    
    prob_matrix = compute_pair_probabilities(sequence)
    
    stems = parse_stems(structure, prob_matrix)
    
    floppy_count = sum(1 for s in stems if s["flag"] == "floppy")
    soft_count = sum(1 for s in stems if s["flag"] == "soft")
    firm_count = sum(1 for s in stems if s["flag"] == "firm")
    
    if floppy_count > 0:
        verdict = "REDESIGN - Contains floppy stems with low ensemble support"
    elif soft_count > len(stems) * 0.5:
        verdict = "NEED PROBING - Multiple stems have moderate ensemble support"
    else:
        verdict = "TRUST - Strong ensemble support for most stems"
    
    heatmap_path = output_dir / "pair_probabilities.png"
    create_heatmap(prob_matrix, heatmap_path)
    
    result = {
        "sequence": sequence,
        "header": header,
        "length": len(sequence),
        "gc_content": gc_content,
        "structure": structure,
        "mfe_energy": mfe_energy,
        "stems": stems,
        "verdict": verdict,
        "metadata": metadata,
        "prob_matrix": prob_matrix,
        "heatmap_path": heatmap_path,
        "stem_counts": {
            "firm": firm_count,
            "soft": soft_count,
            "floppy": floppy_count
        }
    }
    
    if generate_html:
        html_path = output_dir / "report.html"
        generate_html_report(result, html_path)
    
    if generate_markdown:
        md_path = output_dir / "report.md"
        generate_markdown_report(result, md_path)
    
    return result
