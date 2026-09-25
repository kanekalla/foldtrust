#!/usr/bin/env python3
"""Run Layer 5 robustness analysis."""

from pathlib import Path
from foldtrust.benchmark.layer5_robustness import run_robustness_analysis
from foldtrust.utils import find_case_directories

def main():
    cases_dir = Path("data/cases")
    output_dir = Path("benchmarks/outputs/layer5_robustness")
    
    case_dirs = find_case_directories(cases_dir)
    
    print(f"Found {len(case_dirs)} cases:")
    for case_dir in case_dirs:
        print(f"  - {case_dir.name}")
    
    run_robustness_analysis(case_dirs, output_dir)
    
    print("\n✓ Layer 5 complete!")
    print(f"Results in: {output_dir}")

if __name__ == "__main__":
    main()
