"""Pytest configuration and fixtures for FoldTrust tests.

Includes autouse fixture to prevent ViennaRNA parameter state leaks between tests.
"""

import pytest


@pytest.fixture(autouse=True)
def reset_vienna_params():
    """
    Reset ViennaRNA to Turner2004 defaults after every test.

    ViennaRNA 2.7.2 has persistent global parameter state that doesn't fully
    reset between calls in the same process. After params_load_RNA_*(), a fresh
    RNA.md() does NOT take effect if the process has already folded with an
    identical md. This causes parameter state to leak between tests.

    This autouse fixture ensures every test starts with clean Turner2004 defaults.
    """
    yield  # Run the test first

    # After test completes, reset to Turner2004
    try:
        import RNA

        RNA.params_load_RNA_Turner2004()
    except ImportError:
        # ViennaRNA not installed, skip
        pass
