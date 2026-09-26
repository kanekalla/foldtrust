# Layer 0: Disease Case Provenance and Verification

## Summary

Layer 0 establishes the scientific foundation for FoldTrust by ensuring that the five disease-relevant RNA windows are correctly sourced from authoritative public sequence databases with precise genomic coordinates. This layer addresses critical provenance issues discovered in the original case definitions.

## Audit Findings

The original case sequences in `data/cases/` contained significant errors:

1. **SMN2 (`smn2-iss-n1`)**: Sequence was invented — does not occur in RefSeqGene NG_008728.1 or RefSeq NM_017411.4
2. **CFTR (`cftr-5utr`)**: Sequence was invented — does not occur in NM_000492.4 or NG_016465.4
3. **MAPT (`mapt-e10`)**: Sequence was invented — does not occur in NG_007398.2 or NM_001123066.4
4. **HCV (`hcv-ires-dii`)**: Real sequence but wrong region — original was domain III (into core CDS at AF009606.1:148-415), not domain II (5' NTR 44-118)
5. **SARS-CoV-2 (`sars2-fse`)**: Correct sequence and coordinates (NC_045512.2:13462-13542)

**Critical mistake in previous attempts**: NCBI RefSeqGene records use different exon numbering than classic literature. For example, NCBI exon 8 in NG_008728.1 is classic SMN2 exon 7; NCBI exon 11 in NG_007398.2 is classic MAPT exon 10.

## Methodology

### Sequence Extraction

1. Fetch GenBank records from NCBI using E-utilities API
2. Extract feature table annotations to identify exon boundaries
3. Map classic numbering to NCBI numbering using "alternate designation" fields
4. Extract specified 1-based inclusive coordinate ranges
5. Convert T→U for RNA representation
6. Verify strand orientation (all five cases are (+) strand)

### Verification

Script `scripts/verify_cases.py`:
- Fetches each accession from NCBI
- Extracts declared coordinate slice
- Compares with case sequence (T/U normalized)
- Validates header and meta.yaml coordinate agreement
- Reports PASS/FAIL per case

Test suite `tests/test_verify_cases.py`:
- Offline test using cached fixtures in `tests/fixtures/layer0/`
- Online test with live NCBI fetches (skipped when `FOLDTRUST_OFFLINE=1`)
- SHA256 checksums for fixture integrity

## Cases

### SMN2 — Spinal Muscular Atrophy (ISS-N1 Antisense Target)

**Coordinates:** NG_008728.1:31999-32152(+), 154 nt

**Region:** SMN2 exon 7 (classic numbering, NCBI exon 8) + intron 7 +1..+100

**Biology:** The SMN1 gene produces full-length Survival Motor Neuron protein. SMN2, a nearly identical paralogue, produces mostly truncated protein due to a C→T SNP at exon 7 position +6 that creates an exonic splicing silencer, causing exon 7 skipping. Loss of both SMN1 alleles plus insufficient SMN2 exon 7 inclusion causes spinal muscular atrophy. Nusinersen (Spinraza), an FDA-approved antisense oligonucleotide, targets the intronic splicing silencer N1 (ISS-N1) at intron 7 +10..+27 to block splicing repressors and restore exon 7 inclusion.

**Landmarks:**

| Feature | Record coords | Window coords | Sequence | Note |
|---------|--------------|---------------|----------|------|
| Exon 7 (classic, NCBI exon 8) | 31999..32052 | 1..54 | | C6T SNP at pos 6 (U in RNA) |
| Exon 7 / Intron 7 5'ss | 32052\|32053 | 54\|55 | | Splice junction |
| ISS-N1 core | 32062..32076 | 64..78 | CCAGCAUUAUGAAAG | Silencer element |
| Nusinersen target | 32062..32079 | 64..81 | | ASO binds +10..+27 |

**References:**
- Singh NK et al. 2006. Mol Cell Biol 26:1333. DOI: [10.1128/MCB.26.4.1333-1346.2006](https://doi.org/10.1128/MCB.26.4.1333-1346.2006) — ISS-N1 identification
- Finkel RS et al. 2017. N Engl J Med 377:1723. DOI: [10.1056/NEJMoa1702752](https://doi.org/10.1056/NEJMoa1702752) — Nusinersen trial

---

### CFTR — Cystic Fibrosis (5' UTR and Start Codon Region)

**Coordinates:** NM_000492.4:1-200(+), 200 nt

**Region:** 5' UTR (1-70) + CDS start (71-200)

**Biology:** Cystic fibrosis is caused by mutations in the CFTR gene. While most research focuses on coding mutations (e.g., ΔF508), the 5' UTR influences translation efficiency. RNA secondary structure around the start codon affects ribosome scanning and AUG recognition. This window spans the complete annotated 5' UTR (70 nt) plus the first 130 nt of coding sequence to capture start-codon accessibility and ribosome loading context.

**Landmarks:**

| Feature | Record coords | Window coords | Sequence | Note |
|---------|--------------|---------------|----------|------|
| 5' UTR | 1..70 | 1..70 | | Complete annotated UTR |
| Upstream stop codon | 5..7 | 5..7 | UAA | In-frame stop |
| Start codon (ATG) | 71..73 | 71..73 | AUG | Translation start |
| CDS (partial) | 71..200 | 71..200 | | First 130 nt |

**References:**
- Riordan JR et al. 1989. Science 245:1066. DOI: [10.1126/science.2475911](https://doi.org/10.1126/science.2475911) — CFTR identification
- Zielenski J et al. 1991. Genomics 10:214. DOI: [10.1016/0888-7543(91)90503-7](https://doi.org/10.1016/0888-7543(91)90503-7) — Gene structure

---

### MAPT — Frontotemporal Dementia (Exon 10 5' Splice Site Stem-Loop)

**Coordinates:** NG_007398.2:120818-121000(+), 183 nt

**Region:** Intron 9 last 30 nt + MAPT exon 10 (classic numbering, NCBI exon 11) + intron 10 +1..+60

**Biology:** The MAPT gene encodes tau protein. Alternative splicing of exon 10 determines 4-repeat (4R, exon 10 included) vs. 3-repeat (3R, exon 10 skipped) tau isoforms. A stem-loop structure at the exon 10 / intron 10 5' splice site sequesters the splice donor, reducing exon 10 inclusion. Mutations that stabilize the stem (e.g., intron +3, +14, +16 mutations) increase exon 10 inclusion, elevating the 4R/3R ratio and causing frontotemporal dementia with parkinsonism linked to chromosome 17 (FTDP-17).

**Landmarks:**

| Feature | Record coords | Window coords | Sequence | Note |
|---------|--------------|---------------|----------|------|
| Exon 10 (classic, NCBI exon 11) | 120848..120940 | 31..123 | GUGCAGAUAAUUAAUAAGAAG... | 93 nt, encodes R2 repeat |
| Exon/intron junction | 120940\|120941 | 123\|124 | | 5' splice site |
| Stem-loop core | 120935..120959 | 118..142 | GGCAGUGGCCGGUGGGGGCAAGGUG | Exon -5..-1 \| intron +1..+19 |
| S305 codon | 120938..120940 | 121..123 | AGU | S305N at exon -2 |
| Intron +3 | 120943 | 126 | | FTDP-17 G→A |
| Intron +14 | 120954 | 137 | | FTDP-17 C→T |
| Intron +16 | 120956 | 139 | | FTDP-17 C→T |

**References:**
- Grover A et al. 1999. J Biol Chem 274:15134. DOI: [10.1074/jbc.274.21.15134](https://doi.org/10.1074/jbc.274.21.15134) — Stem-loop regulation
- Varani L et al. 1999. Proc Natl Acad Sci USA 96:8229. DOI: [10.1073/pnas.96.14.8229](https://doi.org/10.1073/pnas.96.14.8229) — NMR structure
- Hutton M et al. 1998. Nature 393:702. DOI: [10.1038/31508](https://doi.org/10.1038/31508) — FTDP-17 mutations

---

### HCV — Hepatitis C (IRES Domain II)

**Coordinates:** AF009606.1:44-118(+), 75 nt

**Region:** HCV 5' NTR domain II (standard numbering nt 44-118)

**Biology:** The hepatitis C virus (HCV, genotype 1a strain H77) uses an internal ribosome entry site (IRES) in its 341-nt 5' non-translated region to initiate cap-independent translation. The IRES comprises four domains (I–IV). Domain II forms a conserved stem-loop structure that contributes to IRES architecture and ribosome recruitment. The window spans domain II entirely, from the basal helix (nt 44-47 pairing with 115-118) through the apical loops, and stays within the 5' NTR (the CDS starts at nt 342).

**Landmarks:**

| Feature | Record coords | Window coords | Sequence | Note |
|---------|--------------|---------------|----------|------|
| Domain II basal helix 5' | 44..47 | 1..4 | CCUG | Base pairs with 3' strand |
| Domain II basal helix 3' | 115..118 | 72..75 | CAGG | Pairs with nt 44-47 |
| Domain II (full) | 44..118 | 1..75 | | Complete structured domain |

**Note:** The 5' NTR spans nt 1-341; the AUG start codon is at 342. This window does not overlap coding sequence.

**References:**
- Honda M et al. 1999. J Virol 73:1165. DOI: [10.1128/JVI.73.2.1165-1174.1999](https://doi.org/10.1128/JVI.73.2.1165-1174.1999) — Domain II = nt 44-118
- Lukavsky PJ et al. 2003. Nat Struct Biol 10:1033. DOI: [10.1038/nsb1004](https://doi.org/10.1038/nsb1004) — NMR structure

---

### SARS-CoV-2 — COVID-19 (Frameshifting Element)

**Coordinates:** NC_045512.2:13462-13542(+), 81 nt

**Region:** ORF1a/ORF1ab frameshifting element (slippery site + 3-stem pseudoknot)

**Biology:** SARS-CoV-2 uses a programmed −1 ribosomal frameshift to control the ratio of ORF1a to ORF1ab polyproteins. The frameshift site contains a heptanucleotide slippery sequence (UUUAAAC, nt 13462-13468) followed by a downstream RNA pseudoknot. The pseudoknot comprises three stems that induce ribosomal pausing and slippage. Efficient frameshifting is essential for viral replication, making the FSE a target for antiviral small molecules and antisense oligonucleotides.

**Landmarks:**

| Feature | Record coords | Window coords | Sequence | Note |
|---------|--------------|---------------|----------|------|
| Slippery sequence | 13462..13468 | 1..7 | UUUAAAC | −1 frameshift site |
| Pseudoknot stem 1 | ~13476..13487 | ~15..26 | | Approximate; see structure |
| Pseudoknot stem 2 | ~13488..13498 | ~27..37 | | |
| Pseudoknot stem 3 | ~13519..13542 | ~58..81 | | |

**References:**
- Kelly JA et al. 2020. J Biol Chem 295:10741. DOI: [10.1074/jbc.AC120.013449](https://doi.org/10.1074/jbc.AC120.013449) — FSE mechanism
- Bhatt PR et al. 2021. Science 372:eabf3546. DOI: [10.1126/science.abf3546](https://doi.org/10.1126/science.abf3546) — High-resolution structure

---

## Old vs. New Coordinates

| Case | Old (invented/wrong) | New (verified) | Change |
|------|---------------------|----------------|--------|
| smn2-iss-n1 | No provenance | NG_008728.1:31999-32152(+) 154 nt | Rebuilt from RefSeqGene |
| cftr-5utr | No provenance | NM_000492.4:1-200(+) 200 nt | Rebuilt from RefSeq mRNA |
| mapt-e10 | No provenance | NG_007398.2:120818-121000(+) 183 nt | Rebuilt from RefSeqGene |
| hcv-ires-dii | AF009606.1:148-415 (domain III, 268 nt) | AF009606.1:44-118(+) 75 nt | Corrected to domain II |
| sars2-fse | NC_045512.2:13462-13542(+) 81 nt | NC_045512.2:13462-13542(+) 81 nt | Unchanged (was correct) |

## Verification Output

```
================================================================================
FoldTrust Case Provenance Verification
================================================================================

✓ PASS     cftr-5utr            PASS
✓ PASS     hcv-ires-dii         PASS
✓ PASS     mapt-e10             PASS
✓ PASS     sars2-fse            PASS
✓ PASS     smn2-iss-n1          PASS

================================================================================
Results: 5/5 passed

All cases verified successfully!
```

**Command:** `python3 scripts/verify_cases.py`

## DOI Verification

All DOIs checked via Crossref API (`scripts/check_case_dois.py`):

```
✓ 10.1128/MCB.26.4.1333-1346.2006 — Singh et al. 2006 (SMN2 ISS-N1)
✓ 10.1056/NEJMoa1702752 — Finkel et al. 2017 (Nusinersen)
✓ 10.1126/science.2475911 — Riordan et al. 1989 (CFTR gene)
✓ 10.1016/0888-7543(91)90503-7 — Zielenski et al. 1991 (CFTR structure)
✓ 10.1074/jbc.274.21.15134 — Grover et al. 1999 (MAPT stem-loop)
✓ 10.1073/pnas.96.14.8229 — Varani et al. 1999 (MAPT NMR)
✓ 10.1038/31508 — Hutton et al. 1998 (FTDP-17)
✓ 10.1128/JVI.73.2.1165-1174.1999 — Honda et al. 1999 (HCV domain II)
✓ 10.1038/nsb1004 — Lukavsky et al. 2003 (HCV IRES NMR)
✓ 10.1074/jbc.AC120.013449 — Kelly et al. 2020 (SARS-CoV-2 FSE)
✓ 10.1126/science.abf3546 — Bhatt et al. 2021 (SARS-CoV-2 FSE structure)
```

All DOIs resolve to the cited papers with matching first authors and publication years.

## Case Metrics (ViennaRNA 2.7.2, Turner 2004, 37°C)

| Case | Coordinates | Length | GC% | MFE (kcal/mol) | Firm | Soft | Floppy | Verdict |
|------|-------------|--------|-----|----------------|------|------|--------|---------|
| smn2-iss-n1 | NG_008728.1:31999-32152(+) | 154 | 30.5 | −31.30 | 2 | 6 | 1 | REDESIGN |
| cftr-5utr | NM_000492.4:1-200(+) | 200 | 52.5 | −60.60 | 3 | 1 | 6 | REDESIGN |
| mapt-e10 | NG_007398.2:120818-121000(+) | 183 | 48.1 | −52.30 | 5 | 2 | 4 | REDESIGN |
| hcv-ires-dii | AF009606.1:44-118(+) | 75 | 54.7 | −23.80 | 2 | 3 | 0 | NEED PROBING |
| sars2-fse | NC_045512.2:13462-13542(+) | 81 | 53.1 | −26.00 | 1 | 2 | 1 | REDESIGN |

**Note:** Stem tiers defined by minimum pair probability within each stem (≥3 consecutive base pairs):
- **FIRM**: ≥ 0.85 (high ensemble support)
- **SOFT**: 0.5–0.85 (moderate support)
- **FLOPPY**: < 0.5 (low support, fluctuates in ensemble)

All MFE values match the expected sanity checks specified in task requirements.

**Command:** `python3 scripts/compute_case_metrics.py`

## Known Downstream Staleness

The following files depend on the old (invented/wrong) sequences and coordinates but are outside Layer 0 scope:

1. **`src/foldtrust/benchmark/layer5_robustness.py`**: Hard-coded `CASE_GENOMIC_COORDS` dictionary contains the old wrong coordinates. All Layer 5 parameter-robustness results were computed on the old sequences.

2. **Layer 5 outputs** in `benchmarks/outputs/layer5*/`: CSV/JSON files and figures showing parameter sensitivity (Turner2004 vs. Andronescu2007 vs. Langdon2018) were computed on the old sequences. These results are now invalid and need regeneration.

3. **`docs/benchmark/layer5_*.md`**: Documentation referencing old MFE values and coordinates.

These files are owned by the Layers 1-3 agent and should not be edited by Layer 0. They are flagged here for future correction.

## Limitations

1. **Offline testing**: The `--from-cache` option in `verify_cases.py` and the offline pytest fixture allow verification without network access, but the cached slices must be updated if coordinates change.

2. **Strand orientation**: All five cases are on the (+) strand of their respective reference sequences. The verification script supports (−) strand via reverse-complement but this is untested on real cases.

3. **Exon numbering**: NCBI RefSeqGene records use internal numbering that differs from classic literature numbering (e.g., SMN2 "exon 8" = classic exon 7). The "alternate designation" field in GenBank feature tables resolves this, but requires careful manual mapping.

4. **DOI checking**: The `check_case_dois.py` script uses simple first-author substring matching. It does not check full author lists, journal names, or page numbers. Manual verification of each DOI link is still recommended.

5. **Biological context**: The landmark tables provide coordinate-level precision but do not capture all biological nuance (e.g., MAPT mutations at intron +13 are also pathogenic but not listed as a landmark because the window ends at intron +60).

## Deliverables

- ✅ Five corrected case files (`data/cases/*/sequence.fa`, `meta.yaml`)
- ✅ Verification script (`scripts/verify_cases.py`)
- ✅ Test suite (`tests/test_verify_cases.py`) with cached fixtures (`tests/fixtures/layer0/`)
- ✅ Regenerated outputs:
  - `examples/out/<case>/{report.md, report.html, pair_probabilities.png}`
  - `benchmarks/outputs/layer0_cases/case_metrics.{csv,json}`
  - `tests/fixtures/vienna_api_expected.json`
- ✅ This documentation (`docs/benchmark/layer0_cases.md`)
