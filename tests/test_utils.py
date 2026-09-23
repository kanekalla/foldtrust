"""Tests for utility functions."""

from pathlib import Path
from foldtrust.utils import compute_gc_content, format_structure_ascii


def test_gc_content():
    """Test GC content calculation."""
    assert compute_gc_content("ATCG") == 50.0
    assert compute_gc_content("AAAA") == 0.0
    assert compute_gc_content("GGCC") == 100.0
    assert compute_gc_content("ATCGATCG") == 50.0
    assert compute_gc_content("") == 0.0


def test_format_structure_ascii():
    """Test ASCII structure formatting."""
    seq = "AAAACCCCGGGGTTTT"
    struct = "((((....))))(())"
    
    result = format_structure_ascii(seq, struct, width=10)
    
    assert "AAAACCCCGG" in result
    assert "GGTTTT" in result
    assert "((((....)" in result
