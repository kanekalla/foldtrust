## Layout

- `archiveII/` — ArchiveII (Mathews lab) original tarball + derived JSONL of all 3975 .ct structures.
- `tool_data_releases/mxfold2_zenodo_4430150/` — MXfold2 datasets (Zenodo), split into 3 parts (<20 MB each). Includes ArchiveII (3966 bpseq), bpRNA-1m TR0/VL0/TS0 (canonical pairs only), bpRNA-new, TrainSetA/TestSetA/B, RNAStrAlign.
- `tool_data_releases/spot_rna/` — SPOT-RNA PDB-derived dataset (TR1/VL1/TS1/TS2).
- `rfam/raw/` — 19 Rfam seed Stockholm alignments; `rfam/rfam_seed_ss.jsonl` = projected per-sequence structures (primary); `rfam/rfam_seed_ss.round_angle_only.jsonl` = literal-instruction variant.
- `shape/` — SARS-CoV-2 SHAPE/icSHAPE reactivity CSVs + README/citations + structure tracks, from DasLab/SARS_CoV-2_shape_comparison pinned at commit 75cb3151e200dd7ea62600db50026ce179181fa3 (committed 2021-06-22 19:53 MT).
- `genome/` — NC_045512.2 FASTA + GenBank, `fse_coordinates.json`.
- `vienna_params/` — ViennaRNA rna_turner2004/andronescu2007/langdon2018 .par files (pinned at ViennaRNA commit 1ffec79f5e258896160f7362ced8263450f371dc) + API check result.
- `scripts/` — all conversion/verification scripts.

## Reference-structure record counts

- ArchiveII (`archiveII/archiveII.jsonl.gz`, 3975 records): 16s 110, 23s 35, 5s 1283, RNaseP 454, grp1 98, grp2 11, srp 928, tRNA 557, telomerase 37, tmRNA 462. (16s/23s counts include both full-length files and their `_domainN` files; the MXfold2 copy of archiveII has 3966 entries.)
- bpRNA-1m (via MXfold2 release): TR0 10814, VL0 1300, TS0 1305 (canonical pairs); bpRNA-new 5401.
- SPOT-RNA PDB set: TR1 120, VL1 30, TS1 67, TS2 39 label files.
- Rfam seed projection (`rfam/rfam_seed_ss.jsonl`): **468 records**, 19 families; 25 per family except RF00234 glmS (18: only 18 seed sequences).

| family | Rfam ID | seed seqs | kept | skipped non-ACGU | skipped >400 nt | mean pairs (primary) |
|---|---|---:|---:|---:|---:|---:|
| RF00005 | tRNA | 954 | 25 | 1 | 0 | 20.5 (variant 20.5) |
| RF00001 | 5S_rRNA | 712 | 25 | 17 | 0 | 30.6 (variant 30.6) |
| RF00017 | Metazoa_SRP | 91 | 25 | 0 | 0 | 93.1 (variant 93.1) |
| RF01854 | Bacteria_large_SRP | 92 | 25 | 0 | 0 | 54.4 (variant 52.2) |
| RF00010 | RNaseP_bact_a | 458 | 25 | 8 | 125 | 96.0 (variant 69.8) |
| RF00003 | U1 | 100 | 25 | 0 | 0 | 41.6 (variant 41.6) |
| RF00004 | U2 | 208 | 25 | 2 | 0 | 39.8 (variant 39.8) |
| RF00015 | U4 | 177 | 25 | 2 | 0 | 24.4 (variant 24.4) |
| RF00026 | U6 | 188 | 25 | 1 | 0 | 4.7 (variant 4.7) |
| RF00050 | FMN | 146 | 25 | 0 | 0 | 25.6 (variant 18.5) |
| RF00162 | SAM | 457 | 25 | 0 | 0 | 29.2 (variant 29.2) |
| RF00167 | Purine | 133 | 25 | 0 | 0 | 19.5 (variant 19.5) |
| RF00059 | TPP | 115 | 25 | 1 | 0 | 24.0 (variant 19.4) |
| RF00174 | Cobalamin | 434 | 25 | 1 | 0 | 38.1 (variant 29.8) |
| RF00504 | Glycine | 44 | 25 | 0 | 0 | 20.0 (variant 20.0) |
| RF00380 | Magnesium | 159 | 25 | 1 | 0 | 41.4 (variant 21.8) |
| RF00168 | Lysine | 47 | 25 | 0 | 0 | 46.7 (variant 46.7) |
| RF00234 | glmS | 18 | 18 | 0 | 0 | 36.2 (variant 36.2) |
| RF00023 | tmRNA | 477 | 25 | 9 | 25 | 65.2 (variant 58.7) |

Rfam conversion rules (scripts/rfam_sto_to_jsonl.py): gap columns (. - _ ~) removed; a consensus pair is kept only if both columns are residues and the pair is AU/UA/GC/CG/GU/UG; T->U; sequences containing non-ACGU letters are skipped; sequences >400 nt skipped; filters applied first, then records sorted by seq_id (string sort) and the first 25 kept. Every record verified: equal sequence/structure length, ACGU only, only `().`, balanced, every pair canonical/wobble.

**WUSS note (deviation, please read):** the task said "strip pseudoknot brackets other than <> and ()". In WUSS, `[]` and `{}` are *not* pseudoknots: they are nested helix brackets (same as `()`/`<>`, used for outer helices); pseudoknots are marked with letter pairs (`Aa`, `Bb`, ...). The primary file therefore keeps `()`, `<>`, `[]`, `{}` as nested pairs (parsed with one type-matched stack, which succeeded for all 19 families, proving they are nested) and strips only the letter pseudoknots. The literal variant `rfam_seed_ss.round_angle_only.jsonl` drops `[]`/`{}` as instructed; this removes real helices in RF01854, RF00010, RF00050, RF00059, RF00174, RF00380, RF00023 (see mean-pair column). Note RF00026 (U6) SS_cons has only 5 consensus pairs, so its records are nearly unstructured — a property of the Rfam annotation, not a bug.

## SHAPE files and citations (per DasLab README)

| file | lab / paper |
|---|---|
| incarnato_invitro_reactivity.csv, incarnato_invivo_reactivity.csv | Manfredonia I, et al. NAR 2020;48:12436-12452 (Incarnato lab; README explicitly labels the in vitro set; the in vivo file is from the same lab/paper) |
| pyle_reactivity.csv | Huston NC, et al. Mol Cell 2021;81:584-598 (Pyle lab, in vivo SHAPE-MaP) |
| zhang_invitro_reactivity.csv, zhang_invivo_reactivity.csv | Sun L, et al. Cell 2021;184:1865-1883 (Zhang lab, icSHAPE) |
| secondary_structure_tracks/rnaz_secondary_structures.csv, sarsr_conservation.csv | Rangan R, et al. RNA 2020;26:937-959 |

Format: each reactivity CSV has 29903 lines, one value per nucleotide of NC_045512.2 (DasLab refseq.txt verified identical to NC_045512.2), no header. Missing = `nan`/`NaN` (Incarnato, Zhang) or `-999` (Pyle). Values are **not** uniformly normalized: Incarnato 0–15.2 / 0–9.0, Pyle −14.0–45.7, Zhang 0–1. descriptions.txt describes the normalization used for the heatmap (outlier removal, top-10% scaling, winsorize 0–2, /2), apparently applied in the notebook, not in these CSVs.

## SARS-CoV-2 frameshift element (NC_045512.2, 1-based inclusive)

- Slippery site UUUAAAC: **13462–13468** (unique occurrence between 13000 and 14000; 9 occurrences genome-wide).
- ORF1ab CDS ribosomal slippage: join(266..13468,13468..21555).
- GenBank stem_loop features (Rfam RF00507 inference): stem-loop 1 13476–13503; stem-loop 2 13488–13542.
- FSE span slippery site → end of stem-loop 2: 13462–13542 (81 nt): UUUAAACGGGUUUGCGGUGUAAGUGCAGCCCGUCUUACACCGUGCGGCACAGGCACUAGUACUGAUGUCGUAUACAGGGCU

## ViennaRNA parameter check (pip ViennaRNA 2.7.2)

The Python API exposes `RNA.params_load_RNA_Turner2004`, `RNA.params_load_RNA_Andronescu2007`, `RNA.params_load_RNA_Langdon2018` (and `RNA.params_load(file)`). Loading them **does** change the MFE — but only if a new `RNA.md()` is created after loading and passed to `RNA.fold_compound(seq, md)`. **Pitfall:** in the same process, `RNA.fold(seq)` or `RNA.fold_compound(seq)` without a fresh md kept using the previously cached parameter set, so a naive check reports identical energies for all three sets (seen during this build).

| sequence | Turner2004 | Andronescu2007 | Langdon2018 |
|---|---:|---:|---:|
| RF00005 AAEU02004214.1/433-498 (66 nt) | −10.20 | −10.70 | −12.00 |
| SARS-CoV-2 FSE 13462–13542 (81 nt) | −26.00 | −22.26 | −24.70 |

(kcal/mol.) Built-in loaders and the downloaded .par files gave identical results for every set. Full structures in vienna_params/param_check_result.json.

## Sources tried: succeeded / failed

Succeeded
- ArchiveII via Wayback Machine raw capture of https://rna.urmc.rochester.edu/pub/archiveII.tar.gz (capture 2022-10-17; SHA1 matches Wayback CDX digest).
- ArchiveII (3966) + bpRNA-1m TR0/VL0/TS0 + bpRNA-new + others via MXfold2 Zenodo record 4430150 (found with Zenodo API query `q=MXfold2`, DOI from MXfold2 docs).
- SPOT-RNA PDB_dataset.zip (Dropbox link in SPOT-RNA README).
- Rfam seed alignments for all 19 families from https://rfam.org/family/<acc>/alignment/stockholm.
- DasLab SHAPE repo files via raw.githubusercontent.com pinned to commit 75cb3151… (git blob SHA1 of every file verified against the GitHub tree).
- NCBI efetch FASTA + GenBank for NC_045512.2.
- ViennaRNA .par files via raw.githubusercontent.com (master and pinned commit 1ffec79f… identical); pip ViennaRNA 2.7.2 installed in /workspace/.venv-vienna (not part of the bundle).

Failed / excluded
- https://rna.urmc.rochester.edu/pub/archiveII.tar.gz (and .zip, .tgz, /pub/, /archiveII.tar.gz, RNAStralign.tar.gz): HTTP 404 (the lab's publications page still links pub/archiveII.tar.gz). http:// variant: 502.
- Zenodo API search `q=ArchiveII`: no record containing ArchiveII directly (only unrelated/huge archives, e.g. 22 GB "RiboGRAM Archives").
- Requested Rfam URL form `https://rfam.org/family/<acc>/alignment?acc=<acc>&format=stockholm&download=1`: HTTP 404 {"detail":"Not found."} (new Rfam site); rfam.xfam.org has a TLS hostname mismatch. `-k` was not needed anywhere; no TLS verification was disabled for any bundled file.
- bpRNA-1m official site (bprna.cgrb.oregonstate.edu): downloadable but too large for the budget — dbnFiles.zip 45,237,373 bytes (SHA256 2f6c81133bb2ab23e6baace949c6b3b3e3500bcb73c5fc3d12a6d171d40da94b, 102,318 files, downloaded and verified readable, NOT bundled), fastaFiles.zip 41 MB, bpseqFiles.zip 151 MB, bpRNA_1m_90.zip 4.0 GB, bpRNA_1m.zip 57 GB. bpRNA-1m TR0/VL0/TS0 are available in the MXfold2 release instead.
- SPOT-RNA bpRNA_dataset.zip (https://www.dropbox.com/s/w3kc4iro8ztbf3m/bpRNA_dataset.zip?dl=1, 11,098,220 bytes; TR0 10814/VL0 1300/TS0 1305 — same counts as MXfold2's bpRNA set, but with non-canonical pairs): downloadable, excluded to stay under 40 MB.
- UFold data: Google Drive folder only (repo `data/` holds just a Readme) — not fetched. E2Efold: README has no direct data link found. LinearPartition: no data link in README. GitHub REST API from the box was rate-limited (unauthenticated); commit SHAs/trees were obtained through the GitHub connector instead.
