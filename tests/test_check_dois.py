"""Test DOI resolution checking script."""

import os
import tempfile
from pathlib import Path

from scripts.check_dois import check_doi_resolution, extract_dois_from_file


def test_doi_extraction_basic():
    """Test basic DOI extraction from markdown."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".md", delete=False) as f:
        f.write("Test paper (Smith et al., Nature 1998; doi:10.1038/31508)\n")
        f.write("Another paper doi:10.1186/1748-7188-6-26\n")
        temp_path = f.name

    try:
        dois = extract_dois_from_file(Path(temp_path))
        assert len(dois) == 2
        assert dois[0][0] == "10.1038/31508"
        assert dois[1][0] == "10.1186/1748-7188-6-26"
    finally:
        os.unlink(temp_path)


def test_doi_extraction_markdown_links():
    """Test DOI extraction from markdown links."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".md", delete=False) as f:
        f.write("Reference: [DOI:10.1038/31508](https://doi.org/10.1038/31508)\n")
        f.write("Link: [Smith et al. 1998](https://doi.org/10.1186/1748-7188-6-26)\n")
        temp_path = f.name

    try:
        dois = extract_dois_from_file(Path(temp_path))
        # Should extract DOIs from both link text and URLs
        doi_list = [d[0] for d in dois]
        assert "10.1038/31508" in doi_list
        assert "10.1186/1748-7188-6-26" in doi_list
    finally:
        os.unlink(temp_path)


def test_doi_extraction_parentheses():
    """Test DOI extraction with parentheses in DOI."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".md", delete=False) as f:
        f.write(
            "Reference: [Zielenski et al. 1991](https://doi.org/10.1016/0888-7543(91)90503-7)\n"
        )
        temp_path = f.name

    try:
        dois = extract_dois_from_file(Path(temp_path))
        doi_list = [d[0] for d in dois]
        assert "10.1016/0888-7543(91)90503-7" in doi_list
    finally:
        os.unlink(temp_path)


def test_crossref_resolution_success():
    """Test successful Crossref resolution."""
    # Use a known-good DOI
    result = check_doi_resolution("10.1038/31508")
    assert result["valid"] is True
    assert result["method"] == "crossref"
    assert "first_author" in result
    assert "title" in result


def test_datacite_resolution_success():
    """Test successful DataCite resolution via doi.org redirect."""
    # Use a known Zenodo DOI
    result = check_doi_resolution("10.5281/zenodo.4430150")
    assert result["valid"] is True
    assert result["method"] == "doi.org-redirect"


def test_doi_resolution_failure():
    """Test DOI resolution failure for non-existent DOI."""
    result = check_doi_resolution("10.9999/invalid.doi.12345")
    assert result["valid"] is False
    assert result["error"] == "doi-not-found"


def test_doi_extraction_multiple_per_line():
    """Test multiple DOIs on same line are all extracted."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".md", delete=False) as f:
        f.write("Multiple: (A et al. doi:10.1038/31508; B et al. doi:10.1186/1748-7188-6-26)\n")
        temp_path = f.name

    try:
        dois = extract_dois_from_file(Path(temp_path))
        assert len(dois) == 2
        doi_list = [d[0] for d in dois]
        assert "10.1038/31508" in doi_list
        assert "10.1186/1748-7188-6-26" in doi_list
        # Both should be on line 1
        assert all(d[2] == 1 for d in dois)
    finally:
        os.unlink(temp_path)


def test_doi_extraction_table_format():
    """Test DOI extraction from table rows."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".md", delete=False) as f:
        f.write("| Case | DOI |\n")
        f.write("|------|-----|\n")
        f.write("| test | doi:10.1038/31508 |\n")
        temp_path = f.name

    try:
        dois = extract_dois_from_file(Path(temp_path))
        assert len(dois) == 1
        assert dois[0][0] == "10.1038/31508"
    finally:
        os.unlink(temp_path)


def test_doi_extraction_list_format():
    """Test DOI extraction from list items."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".md", delete=False) as f:
        f.write("- Smith AB et al. 1998. Nature 393:702. DOI: 10.1038/31508\n")
        f.write("- Jones CD. 2011. Science. DOI: 10.1186/1748-7188-6-26\n")
        temp_path = f.name

    try:
        dois = extract_dois_from_file(Path(temp_path))
        assert len(dois) == 2
        doi_list = [d[0] for d in dois]
        assert "10.1038/31508" in doi_list
        assert "10.1186/1748-7188-6-26" in doi_list
    finally:
        os.unlink(temp_path)
