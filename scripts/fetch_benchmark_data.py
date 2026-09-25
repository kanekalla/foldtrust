#!/usr/bin/env python3
"""
Fetch benchmark data from original URLs and verify SHA256 checksums.

Downloads files listed in data/_cache/MANIFEST.md into data/_cache/,
verifying each file against data/_cache/SHA256SUMS.

Usage:
    python scripts/fetch_benchmark_data.py [--verify-only]

Options:
    --verify-only  Only verify existing files, do not download
"""

import argparse
import hashlib
import sys
from pathlib import Path
from urllib.request import urlretrieve
from urllib.error import URLError


def sha256_file(path: Path) -> str:
    """Compute SHA256 hash of a file."""
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        while chunk := f.read(8192):
            h.update(chunk)
    return h.hexdigest()


def load_checksums(checksums_path: Path) -> dict:
    """Load SHA256SUMS file."""
    checksums = {}
    with open(checksums_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split(maxsplit=1)
            if len(parts) == 2:
                sha256, filepath = parts
                checksums[filepath] = sha256
    return checksums


def parse_manifest(manifest_path: Path) -> list:
    """
    Parse MANIFEST.md to extract file paths and source URLs.
    
    Returns list of (relative_path, url, sha256) tuples.
    """
    files = []
    
    with open(manifest_path) as f:
        in_table = False
        for line in f:
            if line.startswith('| file |'):
                in_table = True
                continue
            if in_table:
                if line.startswith('|---'):
                    continue
                if not line.startswith('| `'):
                    in_table = False
                    continue
                
                # Parse table row: | `path` | size | sha256 | url | citation | notes |
                parts = [p.strip() for p in line.split('|')]
                if len(parts) < 5:
                    continue
                
                filepath = parts[1].strip('`')
                sha256 = parts[3]
                url = parts[4]
                
                # Extract URL from markdown link or plain text
                if url.startswith('http'):
                    # Plain URL
                    url_clean = url.split()[0]
                elif 'http' in url:
                    # Extract from various formats
                    import re
                    match = re.search(r'https?://[^\s\)]+', url)
                    if match:
                        url_clean = match.group(0)
                    else:
                        print(f"Warning: Could not parse URL from: {url}")
                        continue
                else:
                    # Derived file, no URL
                    continue
                
                files.append((filepath, url_clean, sha256))
    
    return files


def download_file(url: str, dest: Path, expected_sha256: str) -> bool:
    """
    Download a file and verify its checksum.
    
    Returns True if successful, False otherwise.
    """
    print(f"Downloading: {dest.name}")
    print(f"  URL: {url}")
    
    try:
        urlretrieve(url, dest)
    except URLError as e:
        print(f"  ERROR: Download failed: {e}")
        return False
    
    # Verify checksum
    actual_sha256 = sha256_file(dest)
    
    if actual_sha256 == expected_sha256:
        print(f"  ✓ SHA256 verified")
        return True
    else:
        print(f"  ERROR: SHA256 mismatch")
        print(f"    Expected: {expected_sha256}")
        print(f"    Got:      {actual_sha256}")
        dest.unlink()
        return False


def verify_file(path: Path, expected_sha256: str) -> bool:
    """Verify an existing file's checksum."""
    if not path.exists():
        return False
    
    actual_sha256 = sha256_file(path)
    return actual_sha256 == expected_sha256


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify-only', action='store_true',
                       help='Only verify existing files, do not download')
    parser.add_argument('--cache-dir', type=Path, default=Path('data/_cache'),
                       help='Cache directory (default: data/_cache)')
    args = parser.parse_args()
    
    cache_dir = args.cache_dir
    manifest_path = cache_dir / 'MANIFEST.md'
    checksums_path = cache_dir / 'SHA256SUMS'
    
    if not manifest_path.exists():
        print(f"ERROR: {manifest_path} not found")
        print("Extract the benchmark data bundle first.")
        sys.exit(1)
    
    if not checksums_path.exists():
        print(f"ERROR: {checksums_path} not found")
        sys.exit(1)
    
    # Load checksums
    checksums = load_checksums(checksums_path)
    print(f"Loaded {len(checksums)} checksums from {checksums_path}")
    
    # Parse manifest
    files = parse_manifest(manifest_path)
    print(f"Found {len(files)} downloadable files in {manifest_path}\n")
    
    if args.verify_only:
        print("=== Verification Mode ===\n")
        
        missing = []
        invalid = []
        valid = []
        
        for filepath, url, expected_sha256 in files:
            full_path = cache_dir / filepath
            
            if not full_path.exists():
                missing.append(filepath)
                print(f"✗ Missing: {filepath}")
            elif verify_file(full_path, expected_sha256):
                valid.append(filepath)
                print(f"✓ Valid: {filepath}")
            else:
                invalid.append(filepath)
                print(f"✗ Invalid checksum: {filepath}")
        
        print(f"\n=== Summary ===")
        print(f"Valid: {len(valid)}/{len(files)}")
        print(f"Missing: {len(missing)}")
        print(f"Invalid: {len(invalid)}")
        
        if missing or invalid:
            sys.exit(1)
        else:
            print("\n✓ All files verified")
            sys.exit(0)
    
    # Download mode
    print("=== Download Mode ===\n")
    
    success = 0
    skipped = 0
    failed = 0
    
    for filepath, url, expected_sha256 in files:
        full_path = cache_dir / filepath
        
        # Check if already exists and valid
        if full_path.exists():
            if verify_file(full_path, expected_sha256):
                print(f"✓ Already exists: {filepath}")
                skipped += 1
                continue
            else:
                print(f"! Invalid checksum, re-downloading: {filepath}")
        
        # Create parent directory
        full_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Download
        if download_file(url, full_path, expected_sha256):
            success += 1
        else:
            failed += 1
    
    print(f"\n=== Summary ===")
    print(f"Downloaded: {success}")
    print(f"Skipped (already valid): {skipped}")
    print(f"Failed: {failed}")
    
    if failed > 0:
        print("\n✗ Some downloads failed")
        sys.exit(1)
    else:
        print("\n✓ All files downloaded and verified")
        sys.exit(0)


if __name__ == '__main__':
    main()
