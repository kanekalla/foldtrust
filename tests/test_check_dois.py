"""Test DOI checking script."""

import json
from pathlib import Path
from unittest.mock import Mock, patch

from scripts.check_dois import check_occurrence


def test_check_dois_detects_wrong_author():
    """Test that wrong author is detected."""
    # Mock Crossref metadata for ViennaRNA paper (10.1186/1748-7188-6-26)
    # Real first author: Lorenz
    metadata = {
        "title": ["ViennaRNA Package 2.0"],
        "author": [{"family": "Lorenz", "given": "Ronny"}],
    }

    # Citation with wrong author
    context = "Smith AB, Jones CD, Brown EF. A completely unrelated title about kinase inhibitors. Cell 2019. doi:10.1186/1748-7188-6-26"

    result = check_occurrence("10.1186/1748-7188-6-26", context, metadata)

    assert not result["valid"], "Should fail on wrong author"
    assert "author-check-failed" in result["error"]


def test_check_dois_detects_wrong_title():
    """Test that wrong title is detected."""
    # Mock Crossref metadata
    metadata = {
        "title": ["Association of missense and 5'-splice-site mutations in tau"],
        "author": [{"family": "Hutton", "given": "M"}],
    }

    # Citation with correct author but wrong title
    context = 'Hutton M, Lendon CL, Rizzu P, Baker M. "Deep learning of protein folding structures in neural networks". Nature 1998. doi:10.1038/31508'

    result = check_occurrence("10.1038/31508", context, metadata)

    assert not result["valid"], "Should fail on wrong title"
    assert "title-check-failed" in result["error"]


def test_check_dois_accepts_correct_citation():
    """Test that correct citation passes."""
    metadata = {
        "title": ["ViennaRNA Package 2.0"],
        "author": [{"family": "Lorenz", "given": "Ronny"}],
    }

    context = "Lorenz R, Bernhart SH, Höner zu Siederdissen C, et al. ViennaRNA Package 2.0. Algorithms Mol Biol 2011. doi:10.1186/1748-7188-6-26"

    result = check_occurrence("10.1186/1748-7188-6-26", context, metadata)

    assert result["valid"], f"Should pass for correct citation: {result}"
    assert result["author_ok"] == "ok"


def test_check_dois_accepts_bare_doi():
    """Test that bare DOI (no author/title context) passes with resolve-only."""
    metadata = {
        "title": ["ViennaRNA Package 2.0"],
        "author": [{"family": "Lorenz", "given": "Ronny"}],
    }

    # Bare DOI in table row
    context = "| tRNA-Phe | doi:10.1186/1748-7188-6-26 | 76 nt |"

    result = check_occurrence("10.1186/1748-7188-6-26", context, metadata)

    assert result["valid"], "Should pass for bare DOI"
    assert result["checked"] == "resolve-only"


def test_check_dois_handles_particle_names():
    """Test handling of particle names like 'van Swieten'."""
    metadata = {
        "title": [
            "Association of missense and 5'-splice-site mutations in tau with the inherited dementia FTDP-17"
        ],
        "author": [{"family": "Hutton", "given": "M"}],
    }

    # Should match "van Swieten" by taking last token
    context = "Hutton M, Lendon CL, Rizzu P, Baker M, van Swieten J, Mann D, Lynch T, Heutink P. Nature 1998. doi:10.1038/31508"

    result = check_occurrence("10.1038/31508", context, metadata)

    assert result["valid"], "Should handle particle names correctly"


def test_check_dois_handles_apostrophes():
    """Test handling of apostrophes in names like D'Souza."""
    metadata = {
        "title": [
            "Missense and silent tau gene mutations cause frontotemporal dementia"
        ],
        "author": [{"family": "D'Souza", "given": "Ian"}],
    }

    # Context with straight apostrophe
    context = "D'Souza I, Poorkaj P, Hong M, Nochlin D. Proc Natl Acad Sci USA 1999. doi:10.1073/pnas.96.10.5598"

    result = check_occurrence("10.1073/pnas.96.10.5598", context, metadata)

    assert result["valid"], "Should handle apostrophes correctly"
