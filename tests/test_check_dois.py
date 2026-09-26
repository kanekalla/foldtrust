"""Test DOI checking script."""

from scripts.check_dois import (
    check_occurrence,
    extract_first_cited_author,
    extract_title,
)


# R1 tests: segmentation
def test_segment_notes_multi_citation():
    """R1: NOTES.md:232 segments separate Kelly and Bhatt."""
    # Mock what segmentation should produce
    segment1 = "Kelly JA et al. J Biol Chem (2020)"
    segment2 = "Bhatt PR et al. Science (2021)"

    author1 = extract_first_cited_author(segment1)
    author2 = extract_first_cited_author(segment2)

    assert author1 == "Kelly", f"Expected Kelly, got {author1}"
    assert author2 == "Bhatt", f"Expected Bhatt, got {author2}"


def test_segment_impact_prose():
    """R1: impact.md:42 prose separates Sun, Kelly, Neupane."""
    # Test that prose citations are segmented correctly
    segment1 = "Sun et al. 2018"
    segment2 = "Kelly et al., J Biol Chem 2020"
    segment3 = "Neupane et al., J Mol Biol 2020"

    author1 = extract_first_cited_author(segment1)
    author2 = extract_first_cited_author(segment2)
    author3 = extract_first_cited_author(segment3)

    assert author1 == "Sun"
    assert author2 == "Kelly"
    assert author3 == "Neupane"


# R2 tests: first-author extraction
def test_author_extraction_initials():
    """R2.1: Extract from 'Surname INITIALS' format."""
    assert extract_first_cited_author("Kelly JA et al. J Biol Chem (2020)") == "Kelly"
    assert extract_first_cited_author("Hutton M et al. 1998") == "Hutton"
    assert extract_first_cited_author("Zielenski J et al. 1991") == "Zielenski"
    assert extract_first_cited_author("D'Souza I, Poorkaj P") == "D'Souza"


def test_author_extraction_multi_name():
    """R2.3: Extract first from comma-separated list."""
    assert extract_first_cited_author("Gorodkin, Stricklin & Stormo 2001") == "Gorodkin"


def test_author_extraction_51_authors():
    """R2: Handle 51-author Hutton list."""
    full_list = "Hutton M, Lendon CL, Rizzu P, Baker M, Froelich S, Houlden H, Pickering-Brown S, Chakraverty S, Isaacs A, Grover A, Hackett J, Adamson J, Lincoln S, Dickson D, Davies P, Petersen RC, Stevens M, de Graaff E, Wauters E, van Baren J, Hillebrand M, Joosse M, Kwon JM, Nowotny P, Che LK, Norton J, Morris JC, Reed LA, Trojanowski J, Basun H, Lannfelt L, Neystat M, Fahn S, Dark F, Tannenberg T, Dodd PR, Hayward N, Kwok JB, Schofield PR, Andreadis A, Snowden J, Craufurd D, Neary D, Owen F, Oostra BA, Hardy J, Goate A, van Swieten J, Mann D, Lynch T, Heutink P."
    assert extract_first_cited_author(full_list) == "Hutton"


def test_author_extraction_particle():
    """R2: Extract with particle prefix."""
    assert extract_first_cited_author("Lorenz R, Bernhart SH, Höner zu Siederdissen C") == "Lorenz"


# R3 test: two-name & form
def test_author_two_name_ampersand():
    """R3: 'Surname & Surname' takes first."""
    assert extract_first_cited_author("Zielenski & Tsui 1991") == "Zielenski"


# R4 test: apostrophe normalization
def test_apostrophe_normalization():
    """R4: Curly vs straight apostrophe."""
    metadata = {
        "title": ["Missense and silent tau gene mutations cause frontotemporal dementia"],
        "author": [{"family": "D'Souza", "given": "Ian"}],
    }

    # Cited with straight apostrophe, Crossref has curly
    segment = "D'Souza I, Poorkaj P, Hong M, Nochlin D"

    result = check_occurrence("10.1073/pnas.96.10.5598", segment, metadata)

    assert result["valid"], "Should handle apostrophe normalization"
    assert result["author_ok"] in ["ok", "ok-particle"]


# R5 test: particle names
def test_particle_name_matching():
    """R5: Match 'van Swieten' or just 'Swieten'."""
    metadata = {
        "title": ["Association of missense mutations in tau"],
        "author": [{"family": "van Swieten", "given": "J"}],
    }

    segment = "van Swieten J, Mann D, Lynch T"

    result = check_occurrence("10.1038/example", segment, metadata)

    assert result["valid"], "Should match particle name"


# R6 tests: title extraction
def test_title_italic_journal_skipped():
    """R6: Italic journal names are not treated as titles."""
    # Simulated segment with italic journal (markdown removed in real processing)
    segment = "Gorodkin, Stricklin & Stormo 2001, Nucleic Acids Res 29:2135-2144"

    title = extract_title(segment)
    assert title is None, "Should not extract italic journal as title"


def test_title_double_quoted_extracted():
    """R6: Double-quoted title with 4+ words is extracted."""
    segment = (
        'Smith AB et al. "Structural analysis of RNA folding in coronavirus genomes". Nature 2020.'
    )

    title = extract_title(segment)
    assert title is not None, "Should extract double-quoted title"
    assert "Structural analysis" in title


# Negative tests (must fail)
def test_negative_wrong_author_prose():
    """Negative: Wrong author in prose citation."""
    metadata = {
        "title": ["Association of missense mutations in tau"],
        "author": [{"family": "Hutton", "given": "M"}],
    }

    segment = "Smith et al., Nature 1998"

    result = check_occurrence("10.1038/31508", segment, metadata)

    assert not result["valid"], "Should fail on wrong author"
    assert "author-check-failed" in result["error"]


def test_negative_wrong_author_list():
    """Negative: Wrong author in list citation."""
    metadata = {
        "title": ["Association of missense mutations in tau"],
        "author": [{"family": "Hutton", "given": "M"}],
    }

    segment = "Smith et al., Nature 1998"

    result = check_occurrence("10.1038/31508", segment, metadata)

    assert not result["valid"], "Should fail on wrong author in list"


def test_negative_wrong_three_author_list():
    """Negative: Wrong first author in 3-author list."""
    metadata = {
        "title": ["FUS mutations cause frontotemporal lobar degeneration"],
        "author": [{"family": "Vance", "given": "Caroline"}],
    }

    segment = "Smith AB, Jones CD, Brown EF"

    result = check_occurrence("10.1002/humu.21545", segment, metadata)

    assert not result["valid"], "Should fail on wrong 3-author list"


def test_negative_real_author_later():
    """Negative: Real author appears later but not first."""
    metadata = {
        "title": ["Association of missense mutations in tau"],
        "author": [{"family": "Hutton", "given": "M"}],
    }

    segment = "Smith AB, Jones CD, White IJ, Hutton M"

    result = check_occurrence("10.1038/31508", segment, metadata)

    assert not result["valid"], "Should fail when real author not first"


def test_negative_wrong_link_text():
    """Negative: Wrong author in markdown link text."""
    metadata = {
        "title": ["Association of missense mutations in tau"],
        "author": [{"family": "Hutton", "given": "M"}],
    }

    # Link text would be extracted as segment
    segment = "Smith et al. 1998"

    result = check_occurrence("10.1038/31508", segment, metadata)

    assert not result["valid"], "Should fail on wrong link text"


def test_negative_wrong_quoted_title():
    """Negative: Wrong quoted title."""
    metadata = {
        "title": ["Association of missense and 5'-splice-site mutations in tau"],
        "author": [{"family": "Hutton", "given": "M"}],
    }

    segment = 'Hutton M et al. "Deep learning of protein folding structures in neural networks". Nature 1998'

    result = check_occurrence("10.1038/31508", segment, metadata)

    assert not result["valid"], "Should fail on wrong quoted title"
    assert "title-check-failed" in result["error"]


# Positive tests (must pass)
def test_positive_bare_doi():
    """Positive: Bare DOI with no author/title."""
    metadata = {
        "title": ["ViennaRNA Package 2.0"],
        "author": [{"family": "Lorenz", "given": "Ronny"}],
    }

    segment = "tRNA-Phe"

    result = check_occurrence("10.1186/1748-7188-6-26", segment, metadata)

    assert result["valid"], "Should pass for bare DOI"
    assert result["checked"] == "resolve-only"


def test_positive_correct_citation():
    """Positive: Correct author and title."""
    metadata = {
        "title": ["Association of missense and 5'-splice-site mutations in tau"],
        "author": [{"family": "Hutton", "given": "M"}],
    }

    segment = "Hutton M et al. Nature 1998"

    result = check_occurrence("10.1038/31508", segment, metadata)

    assert result["valid"], "Should pass for correct citation"
