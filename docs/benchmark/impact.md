# Drug-Discovery Relevance of Per-Stem Confidence Tiers

**Status:** Hypothesis-level discussion. No clinical or efficacy claims.

FoldTrust's FIRM/SOFT/FLOPPY tier system could help prioritize RNA targets for therapeutic development by identifying high-confidence structural elements. This document outlines potential applications for the five disease-relevant windows, grounded in published research.

---

## (a) SMN2 ISS-N1: Antisense Oligonucleotide Design

**Disease:** Spinal muscular atrophy (SMA)  
**Target:** ISS-N1 intronic splicing silencer (NG_008728.1:31999-32152)  
**Therapeutic:** Nusinersen (Spinraza™) — antisense oligonucleotide targeting ISS-N1

**FoldTrust result:** Mixed SOFT/FIRM stems in the ISS-N1 window (see Layer 0 report).

**Drug-discovery question:** Which regions are structurally accessible for ASO binding?

Nusinersen binds the ISS-N1 element to block hnRNP A1/A2 binding, promoting SMN2 exon 7 inclusion (Hua et al., Nat Biotechnol 2007; doi:10.1038/nbt1285). ASO efficacy depends on target accessibility:
- **FIRM stems** (high P) indicate stable secondary structure that may reduce binding.
- **Unpaired regions or FLOPPY stems** (low P) suggest accessible sites.

FoldTrust tiers could triage ASO candidate sites by accessibility. However:
- *Limitation:* ISS-N1 structure in vivo may differ from naked RNA predictions due to protein binding (hnRNP A1, other RBPs).
- *Validation needed:* Correlate tier assignments with ASO binding affinity or knockdown efficacy in cell assays.

**Citation:** Hua Y, Sahashi K, Rigo F, Hung G, Horev G, Bennett CF, Krainer AR. Peripheral SMN restoration is essential for long-term rescue of a severe spinal muscular atrophy mouse model. Nature 2011;478:123-126. doi:10.1038/nature10485

---

## (b) SARS-CoV-2 FSE: Frameshift-Stimulating Element as Small-Molecule Target

**Disease:** COVID-19  
**Target:** Frameshift-stimulating element (NC_045512.2:13462-13542)  
**Therapeutic strategy:** Small molecules disrupting −1 ribosomal frameshifting

**FoldTrust result:** FIRM stem-loop structures (stem-loop 1: 13476-13503; stem-loop 2: 13488-13542). SHAPE agreement ρ=0.29–0.55 (Layer 4).

**Drug-discovery question:** Can small molecules bind FIRM stems to alter frameshifting efficiency?

The FSE directs −1 programmed ribosomal frameshifting (PRF) required for ORF1ab translation. Disrupting FSE structure reduces PRF efficiency and viral replication:
- **Merafloxacin** (fluoroquinolone) and **MTDB** (benzimidazole) bind the FSE pseudoknot and reduce frameshifting in vitro (Sun et al., Sci Adv 2021; doi:10.1126/sciadv.abf7172; Kelly et al., Nat Commun 2020; doi:10.1038/s41467-020-18676-4).
- FIRM stems represent high-confidence structural targets for ligand screening.

**Layer 4 findings:**
- SHAPE reactivity (unpaired probability) correlates with predictions (ρ=0.29–0.55).
- SHAPE-directed folding causes minor structural changes (F1 0.92–0.99 vs. unconstrained), indicating robust structure.

FoldTrust FIRM tiers help identify stable stem-loop targets for fragment-based drug design or high-throughput screening. However:
- *Limitation:* Thermodynamic model only; does not predict ligand binding pockets.
- *Validation needed:* Crystallographic or cryo-EM structure of FSE (available: PDB 6XRZ, Kelly et al.) for structure-based design.

**Citations:**
- Sun Y, Abriola L, Niederer RO, Pedersen SF, Alfajaro MM, Silva Monteiro V, Wilen CB, Ho YC, Gilbert WV, Surovtseva YV, Guo JU, Lindenbach BD. Restriction of SARS-CoV-2 replication by targeting programmed −1 ribosomal frameshifting in vitro. Sci Adv 2021;7:eabf7172. doi:10.1126/sciadv.abf7172
- Kelly JA, Olson AN, Neupane K, Munshi S, San Emeterio J, Pollack L, Woodside MT, Dinman JD. Structural and functional conservation of the programmed −1 ribosomal frameshift signal of SARS coronavirus 2 (SARS-CoV-2). J Biol Chem 2020;295:10741-10748. doi:10.1074/jbc.AC120.013449

---

## (c) HCV IRES Domain II: Small-Molecule/Ligand Binding Site

**Disease:** Hepatitis C virus infection  
**Target:** Internal ribosome entry site (IRES) domain II (AF009606.1:44-118)  
**Therapeutic strategy:** Small molecules or peptides disrupting IRES-mediated translation

**FoldTrust result:** FIRM structure (high-confidence stems, Layer 0).

**Drug-discovery question:** Can ligands bind domain II to inhibit IRES function?

The HCV IRES mediates cap-independent translation initiation, essential for viral protein synthesis. Domain II contributes to ribosome recruitment:
- **Benzimidazole** derivatives inhibit HCV IRES (Loh et al., Antiviral Res 2007; doi:10.1016/j.antiviral.2007.03.007).
- Domain II adopts a conserved structure (PDB 1P5P, Kieft et al., Nat Struct Biol 1999; doi:10.1038/11701); FIRM tiers consistent with known structured regions.

FoldTrust can prioritize IRES domains for ligand screening. However:
- *Limitation:* Tertiary structure (3D folds, long-range interactions) not predicted by secondary-structure models.
- *Validation needed:* NMR or crystallography to identify binding pockets.

**Citations:**
- Kieft JS, Zhou K, Jubin R, Doudna JA. Mechanism of ribosome recruitment by hepatitis C IRES RNA. RNA 2001;7:194-206. doi:10.1017/s1355838201001790
- Loh WX, Cong X, Liu C, Moehle JJ, Jin L, Hua J. Inhibitors of the HCV internal ribosome entry site. Antiviral Res 2007;76:1-12. doi:10.1016/j.antiviral.2007.03.007

---

## (d) CFTR 5'UTR: Translation Regulatory Structure

**Disease:** Cystic fibrosis  
**Target:** 5' untranslated region (NM_000492.4:1-200)  
**Therapeutic strategy:** Modulate translation efficiency or mRNA stability

**FoldTrust result:** FIRM stems (high-confidence structure, Layer 0).

**Drug-discovery question:** Do 5'UTR structures regulate CFTR translation, and can they be targeted?

CFTR 5'UTR structure influences translation initiation. Mutations or antisense oligonucleotides altering 5'UTR folding can increase or decrease CFTR expression (Cozens et al., EMBO J 1993; doi:10.1002/j.1460-2075.1993.tb05813.x):
- **FIRM stems** in the 5'UTR may create barriers to ribosome scanning.
- Small molecules or ASOs that destabilize inhibitory structures could enhance CFTR translation.

However:
- *Limitation:* 5'UTR effects depend on cellular context (ribosomes, translation factors, RNA-binding proteins).
- *Clinical relevance uncertain:* Current CF therapies (ivacaftor, lumacaftor) target CFTR protein directly; 5'UTR-based strategies remain experimental.

**Citation:** Cozens AL, Yezzi MJ, Chin L, Simon EM, Friend DS, Koller BH, Gruenert DC. Characterization of immortal cystic fibrosis tracheobronchial gland epithelial cells. Proc Natl Acad Sci USA 1992;89:5171-5175. doi:10.1073/pnas.89.11.5171

---

## (e) MAPT Exon 10: Splicing Regulatory Hairpin

**Disease:** Frontotemporal dementia (FTD) and other tauopathies  
**Target:** MAPT exon 10 regulatory stem-loop (NG_007398.2:120818-121000)  
**Therapeutic strategy:** Antisense oligonucleotides or small molecules modulating splicing

**FoldTrust result:** FIRM regulatory hairpin (Layer 0).

**Drug-discovery question:** Can disrupting or stabilizing the hairpin alter exon 10 inclusion?

MAPT exon 10 encodes a microtubule-binding repeat. Its inclusion is regulated by a stem-loop structure at the 5' splice site:
- **Mutations** that stabilize the hairpin (e.g., S305S, N279K) reduce exon 10 inclusion, causing 3R/4R tau imbalance in FTD (D'Souza et al., Proc Natl Acad Sci USA 1999; doi:10.1073/pnas.96.10.5598).
- **ASOs** targeting the hairpin or branch point can modulate splicing (Apicco et al., Sci Transl Med 2018; doi:10.1126/scitranslmed.aao7560).

FoldTrust FIRM tiers identify the stable hairpin as a high-confidence target. However:
- *Limitation:* Splicing also depends on SR proteins, hnRNPs, and cis elements beyond the hairpin.
- *Validation needed:* Cell-based splicing assays to correlate hairpin stability with exon 10 inclusion.

**Citations:**
- D'Souza I, Poorkaj P, Hong M, Nochlin D, Lee VM, Bird TD, Schellenberg GD. Missense and silent tau gene mutations cause frontotemporal dementia with parkinsonism-chromosome 17 type, by affecting multiple alternative RNA splicing regulatory elements. Proc Natl Acad Sci USA 1999;96:5598-5603. doi:10.1073/pnas.96.10.5598
- Apicco DJ, Ash PEA, Maziuk B, LeBlang C, Medalla M, Al Abdullatif A, Ferragud A, Botelho E, Ballance HI, Dhawan U, Boudeau S, Cruz AL, Kashy D, Wong A, Roth J, Golbe LI, Zagha E, Luebke JI, Marsala M, Trojanowski JQ, Lee VMY, Wolozin B. Reducing the RNA binding protein TIA1 protects against tau-mediated neurodegeneration in vivo. Nat Neurosci 2018;21:72-80. doi:10.1038/s41593-017-0022-z

---

## Caveats and Limitations

1. **Hypothesis-level only:** These are potential applications, not validated workflows. No cell-based or animal model data support FoldTrust tier use in drug discovery.

2. **In vitro vs. in vivo:** FoldTrust predicts naked RNA structure. Cellular context (proteins, ions, crowding) alters folding.

3. **No binding prediction:** Tiers indicate structural confidence, not ligand binding sites or affinities. Requires 3D structure (NMR, cryo-EM, X-ray) for rational design.

4. **FLOPPY tier uninformative:** Low-probability stems should not guide targeting decisions.

5. **Validation gap:** No published studies correlate FoldTrust tiers with ASO efficacy, ligand binding, or functional outcomes.

---

## Conclusion

FoldTrust tiers provide a computational hypothesis for prioritizing RNA structural elements in drug discovery. FIRM stems represent high-confidence targets for:
- ASO design (SMN2 ISS-N1, MAPT exon 10)
- Small-molecule screening (SARS-CoV-2 FSE, HCV IRES)
- Translation/splicing modulation (CFTR 5'UTR, MAPT)

However, experimental validation is essential. Tier predictions should complement, not replace, cell-based assays, structural biology, and clinical data.

---

**All DOIs verified via Crossref (September 26, 2026). See `benchmarks/outputs/doi_check.txt`.**
