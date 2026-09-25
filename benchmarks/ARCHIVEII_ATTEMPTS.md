"""
ArchiveII Dataset Acquisition Attempts - Documentation

Attempts made to obtain the real ArchiveII dataset:

1. Mathews Lab (rna.urmc.rochester.edu):
   - Direct URL: https://rna.urmc.rochester.edu/RNAstructure/Supplemental/ArchiveII/archiveII.tar.gz
   - Result: HTTP 404 Not Found
   - Note: Main site accessible, but archiveII.tar.gz endpoint does not exist

2. RNA STRAND v2.0:
   - URL: http://www.rnasoft.ca/strand/download/RNAstrand_v2.0.tar.gz
   - Result: Download failed (timeout/connection refused)
   
3. Published benchmark repositories checked:
   - mxfold2 (https://github.com/mxfold/mxfold2): Cloned, no .ct/.bpseq files found
   - LinearPartition (https://github.com/LinearFold/LinearPartition): No structure files
   - allegro (ucrbioinfo): ArchiveII.tar.gz not found
   - E2Efold (e2efold-productive): ArchiveII_all.pkl not found
   - DasLab biers: setG.json not accessible
   
4. Public databases:
   - Rfam: Available but requires family-by-family extraction and parsing
   - PDB: Individual structures available but no pre-compiled ArchiveII equivalent

Conclusion:
The original ArchiveII dataset from the Mathews lab appears to no longer be hosted
at its documented URLs. Published RNA structure prediction papers from 2018-2024
reference ArchiveII but don't redistribute it in their repos (likely due to size).

Without access to the validated ArchiveII dataset, reporting reference structure
accuracy metrics would require either:
a) Using an invented/hand-curated set (not acceptable per instructions)
b) Building a validated set from PDB/Rfam (substantial data engineering beyond scope)

Per user instruction: "If you truly cannot get real reference data, say so and 
do not report reference metrics at all."

Decision: Remove reference structure accuracy from benchmark suite. Focus on:
- Calibration analysis using disease case windows (sars2-fse, etc.)
- SHAPE probing correlation (obtain real data)
- Robustness analysis (parameter sets, temperature, window jitter)
"""