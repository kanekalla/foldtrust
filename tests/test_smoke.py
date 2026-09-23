"""Smoke tests for CLI and end-to-end workflow."""

import pytest
from pathlib import Path
from foldtrust.vienna import check_viennarna


@pytest.mark.skipif(not check_viennarna(), reason="ViennaRNA not installed")
def test_smoke_hairpin(tmp_path):
    """Smoke test with a simple hairpin sequence."""
    from foldtrust.core import process_sequence
    
    hairpin_fasta = tmp_path / "hairpin.fa"
    hairpin_fasta.write_text(">test_hairpin\nGGGAAACCC\n")
    
    output_dir = tmp_path / "output"
    
    result = process_sequence(hairpin_fasta, output_dir)
    
    assert result["length"] == 9
    assert result["structure"] is not None
    assert len(result["stems"]) >= 0
    assert (output_dir / "report.html").exists()
    assert (output_dir / "report.md").exists()
    assert (output_dir / "pair_probabilities.png").exists()


def test_case_structure():
    """Test that all disease cases have required files."""
    cases_dir = Path("data/cases")
    
    if not cases_dir.exists():
        pytest.skip("data/cases not found (run from repo root)")
    
    expected_cases = ["sars2-fse", "smn2-iss-n1", "cftr-5utr", "mapt-e10", "hcv-ires-dii"]
    
    for case_name in expected_cases:
        case_dir = cases_dir / case_name
        assert case_dir.exists(), f"Case {case_name} not found"
        assert (case_dir / "sequence.fa").exists(), f"sequence.fa missing for {case_name}"
        assert (case_dir / "meta.yaml").exists(), f"meta.yaml missing for {case_name}"
