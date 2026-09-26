"""Test DOI checking script."""

from scripts.check_dois import check_occurrence


def test_check_dois_detects_wrong_author_et_al():
    """Test that wrong author in 'et al.' citation is detected."""
    metadata = {
        "title": ["Association of missense and 5'-splice-site mutations in tau"],
        "author": [{"family": "Hutton", "given": "M"}],
    }

    context = "Smith et al., Nature 1998; doi:10.1038/31508"

    result = check_occurrence("10.1038/31508", context, metadata)

    assert not result["valid"], "Should fail on wrong 'et al.' first author"
    assert "author-check-failed" in result["error"]


def test_check_dois_detects_wrong_first_author_three_author_list():
    """Test that wrong first author in 3-author list is detected."""
    metadata = {
        "title": ["FUS mutations cause frontotemporal lobar degeneration"],
        "author": [{"family": "Vance", "given": "Caroline"}],
    }

    context = "Smith AB, Jones CD, Brown EF. FUS mutations cause frontotemporal lobar degeneration. Science 2009. doi:10.1002/humu.21545"

    result = check_occurrence("10.1002/humu.21545", context, metadata)

    assert not result["valid"], "Should fail on wrong first author in 3-author list"
    assert "author-check-failed" in result["error"]


def test_check_dois_detects_real_author_later_in_list():
    """Test that real first author appearing later (but not first) fails."""
    metadata = {
        "title": ["Association of missense and 5'-splice-site mutations in tau"],
        "author": [{"family": "Hutton", "given": "M"}],
    }

    context = "Smith AB, Jones CD, Hutton M, Baker M. Nature 1998. doi:10.1038/31508"

    result = check_occurrence("10.1038/31508", context, metadata)

    assert not result["valid"], "Should fail when real author appears but not first"
    assert "author-check-failed" in result["error"]


def test_check_dois_handles_apostrophes_curly_vs_straight():
    """Test that curly apostrophe (U+2019) in Crossref matches straight in citation."""
    metadata = {
        "title": ["Missense and silent tau gene mutations cause frontotemporal dementia"],
        "author": [{"family": "D'Souza", "given": "Ian"}],
    }

    context = "D'Souza I, Poorkaj P, Hong M, Nochlin D. Proc Natl Acad Sci USA 1999. doi:10.1073/pnas.96.10.5598"

    result = check_occurrence("10.1073/pnas.96.10.5598", context, metadata)

    assert result["valid"], "Should handle curly vs straight apostrophe"
    assert result["author_ok"] in ["ok", "ok-particle"]


def test_check_dois_ignores_prime_symbols_not_as_title_quotes():
    """Test that 5' in prose like '5'UTR' is not treated as title quotes."""
    metadata = {
        "title": ["Association of missense and 5'-splice-site mutations in tau"],
        "author": [{"family": "Hutton", "given": "M"}],
    }

    context = "Hutton M et al. described 5'UTR variants and 5' splice site mutations. Nature 1998. doi:10.1038/31508"

    result = check_occurrence("10.1038/31508", context, metadata)

    assert result[
        "valid"
    ], "Should pass (5' is not treated as a quote delimiter for title extraction)"


def test_check_dois_isolates_neighbouring_entries():
    """Test that neighbouring entry's author list doesn't contaminate context."""
    metadata = {
        "title": ["Structural basis of RNA folding and recognition in coronavirus"],
        "author": [{"family": "Zhang", "given": "Kaiming"}],
    }

    # Zhang et al. followed by Sun et al. in next line
    context = "Zhang K, Zheludev IN, Hagey RJ, Haslecker R, Hou YJ. doi:10.1038/s41594-021-00653-y"

    result = check_occurrence("10.1038/s41594-021-00653-y", context, metadata)

    assert result["valid"], "Should isolate citation from neighbouring entry"
    assert result["author_ok"] in ["ok", "ok-particle"]


def test_check_dois_handles_long_51_author_list():
    """Test that 51-author list (Hutton 1998) is parsed correctly."""
    metadata = {
        "title": ["Association of missense and 5'-splice-site mutations in tau"],
        "author": [{"family": "Hutton", "given": "M"}],
    }

    full_list = "Hutton M, Lendon CL, Rizzu P, Baker M, Froelich S, Houlden H, Pickering-Brown S, Chakraverty S, Isaacs A, Grover A, Hackett J, Adamson J, Lincoln S, Dickson D, Davies P, Petersen RC, Stevens M, de Graaff E, Wauters E, van Baren J, Hillebrand M, Joosse M, Kwon JM, Nowotny P, Che LK, Norton J, Morris JC, Reed LA, Trojanowski J, Basun H, Lannfelt L, Neystat M, Fahn S, Dark F, Tannenberg T, Dodd PR, Hayward N, Kwok JB, Schofield PR, Andreadis A, Snowden J, Craufurd D, Neary D, Owen F, Oostra BA, Hardy J, Goate A, van Swieten J, Mann D, Lynch T, Heutink P."
    context = f"{full_list} Nature 1998. doi:10.1038/31508"

    result = check_occurrence("10.1038/31508", context, metadata)

    assert result["valid"], "Should handle 51-author list without window truncation"
    assert result["author_ok"] in ["ok", "ok-particle"]


def test_check_dois_skips_unquoted_title():
    """Test that unquoted titles are not checked (known loophole per J1 spec)."""
    metadata = {
        "title": ["Association of missense and 5'-splice-site mutations in tau"],
        "author": [{"family": "Hutton", "given": "M"}],
    }

    context = "Hutton M, Lendon CL, Rizzu P. Deep learning for protein structure prediction at genome scale. Nature 1998. doi:10.1038/31508"

    result = check_occurrence("10.1038/31508", context, metadata)

    assert result["valid"], "Should pass - unquoted titles are not checked (loophole)"
    assert result["title_ok"] == "skip-no-title"


def test_check_dois_accepts_correct_citation():
    """Test that correct citation passes."""
    metadata = {
        "title": ["ViennaRNA Package 2.0"],
        "author": [{"family": "Lorenz", "given": "Ronny"}],
    }

    context = "Lorenz R, Bernhart SH, Höner C, Tafer H, Flamm C, Stadler PF. ViennaRNA Package 2.0. Algorithms Mol Biol 2011. doi:10.1186/1748-7188-6-26"

    result = check_occurrence("10.1186/1748-7188-6-26", context, metadata)

    assert result["valid"], f"Should pass for correct citation: {result}"
    assert result["author_ok"] in ["ok", "ok-particle"]


def test_check_dois_accepts_bare_doi():
    """Test that bare DOI (no author/title context) passes with resolve-only."""
    metadata = {
        "title": ["ViennaRNA Package 2.0"],
        "author": [{"family": "Lorenz", "given": "Ronny"}],
    }

    context = "| tRNA-Phe | doi:10.1186/1748-7188-6-26 | 76 nt |"

    result = check_occurrence("10.1186/1748-7188-6-26", context, metadata)

    assert result["valid"], "Should pass for bare DOI"
    assert result["checked"] == "resolve-only"


def test_check_dois_detects_wrong_title():
    """Test that wrong title is detected."""
    metadata = {
        "title": ["Association of missense and 5'-splice-site mutations in tau"],
        "author": [{"family": "Hutton", "given": "M"}],
    }

    context = 'Hutton M, Lendon CL, Rizzu P, Baker M. "Deep learning of protein folding structures in neural networks". Nature 1998. doi:10.1038/31508'

    result = check_occurrence("10.1038/31508", context, metadata)

    assert not result["valid"], "Should fail on wrong title"
    assert "title-check-failed" in result["error"]
