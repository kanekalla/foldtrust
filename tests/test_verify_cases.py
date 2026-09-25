"""Test case provenance verification.

This test suite verifies that:
1. The cached fixtures match their SHA256 checksums
2. Case sequences match the cached fixtures (offline test)
3. Case sequences match live NCBI fetches (online test, skipped when offline)
"""

import hashlib
import os
import sys
from pathlib import Path

import pytest

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

# Import the verify module
verify_module_path = Path(__file__).parent.parent / "scripts" / "verify_cases.py"
spec = __import__("importlib.util").util.spec_from_file_location("verify_cases", verify_module_path)
verify_cases = __import__("importlib.util").util.module_from_spec(spec)
spec.loader.exec_module(verify_cases)


@pytest.fixture
def fixtures_dir():
    """Return path to layer0 fixtures directory."""
    return Path(__file__).parent / "fixtures" / "layer0"


@pytest.fixture
def cases_dir():
    """Return path to cases directory."""
    return Path(__file__).parent.parent / "data" / "cases"


def test_fixture_checksums(fixtures_dir):
    """Verify that cached fixtures match their SHA256 checksums."""
    sums_file = fixtures_dir / "SHA256SUMS"
    assert sums_file.exists(), "SHA256SUMS file not found"
    
    lines = sums_file.read_text().strip().split('\n')
    
    for line in lines:
        expected_hash, filename = line.split(None, 1)
        fixture_file = fixtures_dir / filename
        
        assert fixture_file.exists(), f"Fixture file not found: {filename}"
        
        content = fixture_file.read_bytes()
        actual_hash = hashlib.sha256(content).hexdigest()
        
        assert actual_hash == expected_hash, f"Checksum mismatch for {filename}"


def test_cases_offline(fixtures_dir, cases_dir):
    """Verify cases against cached fixtures (offline test)."""
    cache_dir = Path("data/_cache/ncbi")
    cache_dir.mkdir(parents=True, exist_ok=True)
    
    case_dirs = sorted([d for d in cases_dir.iterdir() if d.is_dir()])
    
    assert len(case_dirs) > 0, "No case directories found"
    
    failed = []
    
    for case_dir in case_dirs:
        case_name = case_dir.name
        passed, message = verify_cases.verify_case(case_dir, cache_dir, from_cache=fixtures_dir)
        
        if not passed:
            failed.append((case_name, message))
    
    if failed:
        fail_msg = "\n".join(f"  {name}: {msg}" for name, msg in failed)
        pytest.fail(f"Cases failed verification:\n{fail_msg}")


def test_cases_online(cases_dir):
    """Verify cases against live NCBI fetches (online test, skipped when offline)."""
    # Skip if explicitly offline
    if os.environ.get('FOLDTRUST_OFFLINE') == '1':
        pytest.skip("Skipping online test: FOLDTRUST_OFFLINE=1")
    
    cache_dir = Path("data/_cache/ncbi")
    cache_dir.mkdir(parents=True, exist_ok=True)
    
    case_dirs = sorted([d for d in cases_dir.iterdir() if d.is_dir()])
    
    assert len(case_dirs) > 0, "No case directories found"
    
    failed = []
    network_error = False
    
    for case_dir in case_dirs:
        case_name = case_dir.name
        
        try:
            passed, message = verify_cases.verify_case(case_dir, cache_dir, from_cache=None)
            
            if not passed:
                # Check if it's a network error
                if "Failed to fetch" in message or "URLError" in message:
                    network_error = True
                else:
                    failed.append((case_name, message))
        except Exception as e:
            error_str = str(e)
            if "URLError" in error_str or "timeout" in error_str.lower():
                network_error = True
            else:
                failed.append((case_name, f"Exception: {e}"))
    
    if network_error:
        pytest.skip("Skipping online test: network error")
    
    if failed:
        fail_msg = "\n".join(f"  {name}: {msg}" for name, msg in failed)
        pytest.fail(f"Cases failed verification:\n{fail_msg}")
