"""
SHAPE Data Investigation for SARS-CoV-2 FSE

Data Source Found:
- Repository: https://github.com/DasLab/SARS_CoV-2_shape_comparison
- Contains: Genome-wide SHAPE reactivity data from multiple studies
  - Zhang in vivo: zhang_invivo_reactivity.csv  (29,903 positions)
  - Zhang in vitro: zhang_invitro_reactivity.csv (29,903 positions)
  - Incarnato in vivo: incarnato_invivo_reactivity.csv (29,902 positions)
  - Incarnato in vitro: incarnato_invitro_reactivity.csv (29,902 positions)
  - Pyle: pyle_reactivity.csv (29,903 positions)

Issue:
The FSE sequence in data/cases/sars2-fse/sequence.fa (181 nt) does NOT match
any substring of the reference genome in the DasLab repository. This suggests:
1. The FoldTrust FSE sequence may be from a different viral isolate/strain
2. The coordinates may need remapping
3. The sequence annotation may be from a different source

To properly implement SHAPE correlation:
1. Need to identify the exact genomic coordinates of the FoldTrust FSE sequence
2. Extract corresponding SHAPE reactivity values from the genome-wide data
3. Map them to the 181 nt window

Alternative approach:
Use one of the disease cases that has published SHAPE data with clear coordinate
mapping, or select a well-documented region of SARS-CoV-2 with known coordinates.

Time constraint decision:
Given the coordinate mismatch and time to properly resolve it, I'll focus on:
1. Implementing proper robustness analysis (parameter sets, temperature, window jitter)
2. Document the SHAPE data availability and coordinate mapping issue
3. Provide infrastructure that could be completed with correct coordinates
"""