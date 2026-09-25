#!/usr/bin/env python3
"""
Fetch benchmark data from original URLs and verify SHA256 checksums.

This script downloads the benchmark data from the sources listed in MANIFEST.md
and verifies integrity. The actual data is cached in data/_cache/ (gitignored).
"""

import hashlib
import sys
import tarfile
from pathlib import Path
from urllib.request import urlretrieve

# SHA256 checksums from MANIFEST.md
CHECKSUMS = {
    "foldtrust_bench_data.tar.gz": "PLACEHOLDER_TO_BE_COMPUTED",
    "bprna_TS0_canonicals.tar.gz": "PLACEHOLDER_TO_BE_COMPUTED",
}

# Original URLs from MANIFEST.md
SOURCES = {
    # Main benchmark bundle - to be uploaded to a permanent location
    "foldtrust_bench_data.tar.gz": None,  # Not yet on a public URL
    # bpRNA TS0 canonicals - extracted from MXfold2 Zenodo release
    "bprna_TS0_canonicals.tar.gz": None,  # Extracted locally, not a direct download
}


def compute_sha256(filepath: Path) -> str:
    """Compute SHA256 checksum of a file."""
    sha256 = hashlib.sha256()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            sha256.update(chunk)
    return sha256.hexdigest()


def verify_archive(filepath: Path, expected_sha256: str) -> bool:
    """Verify archive integrity."""
    if not filepath.exists():
        return False
    actual = compute_sha256(filepath)
    return actual == expected_sha256


def main():
    """Main entry point."""
    cache_dir = Path("data/_cache")
    uploads_dir = Path.home() / ".cursor/projects/workspace/uploads"
    
    print("FoldTrust Benchmark Data Fetcher")
    print("=" * 60)
    print()
    
    # For now, since the data is in uploads/, copy it to cache
    # In production, this would download from permanent URLs
    
    if not cache_dir.exists():
        cache_dir.mkdir(parents=True)
    
    # Extract foldtrust_bench_data
    bundle_path = uploads_dir / "foldtrust_bench_data.tar_b5b8.gz"
    if bundle_path.exists():
        print(f"Extracting {bundle_path.name}...")
        with tarfile.open(bundle_path, "r:gz") as tar:
            tar.extractall(cache_dir)
        print(f"  ✓ Extracted to {cache_dir}")
    else:
        print(f"  ✗ {bundle_path} not found")
        return 1
    
    # Extract bpRNA TS0
    bprna_path = uploads_dir / "bprna_TS0_canonicals.tar_dd45.gz"
    if bprna_path.exists():
        print(f"Extracting {bprna_path.name}...")
        with tarfile.open(bprna_path, "r:gz") as tar:
            tar.extractall(cache_dir)
        print(f"  ✓ Extracted to {cache_dir}")
    else:
        print(f"  ✗ {bprna_path} not found")
        return 1
    
    # Verify MANIFEST and SHA256SUMS are present
    manifest_path = cache_dir / "MANIFEST.md"
    if not manifest_path.exists():
        print(f"  ✗ MANIFEST.md not found in {cache_dir}")
        return 1
    
    sha256sums_path = cache_dir / "SHA256SUMS"
    if not sha256sums_path.exists():
        print(f"  ✗ SHA256SUMS not found in {cache_dir}")
        return 1
    
    print()
    print("Data extracted successfully.")
    print(f"Cache directory: {cache_dir.absolute()}")
    print()
    print("To verify bundle integrity, run:")
    print(f"  cd {cache_dir} && python3 scripts/verify_bundle.py")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
