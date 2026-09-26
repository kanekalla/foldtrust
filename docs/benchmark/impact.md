# Drug-Discovery Relevance of Per-Stem Confidence Tiers

**Status:** Hypothesis-level discussion. No clinical or efficacy claims.

FoldTrust's FIRM/SOFT/FLOPPY tier system could help prioritize RNA targets for therapeutic development by identifying high-confidence structural elements. This document outlines potential applications for the five disease-relevant windows, grounded in published research.

---

## (a) SMN2 ISS-N1: Antisense Oligonucleotide Design

**Disease:** Spinal muscular atrophy (SMA)  
**Target:** ISS-N1 intronic splicing silencer (NG_008728.1:31999-32152)  
**Therapeutic:** Nusinersen (Spinraza™) — antisense oligonucleotide targeting ISS-N1

**FoldTrust result:** SMN2 window: 2 FIRM / 6 SOFT / 1 FLOPPY stems (min-prob tiers; benchmarks/outputs/layer0_cases/case_metrics.csv); Layer 5: FIRM retention 0.40 with 25–100 nt genomic flanks, 0.00 under Langdon2018 (benchmarks/outputs/layer5/*_stem_retention.csv).

**Drug-discovery question:** Which regions are structurally accessible for ASO binding?

Nusinersen binds the ISS-N1 element to block hnRNP A1/A2 binding, promoting SMN2 exon 7 inclusion (Hua et al., Am J Hum Genet 2008; doi:10.1016/j.ajhg.2008.01.014; ISS-N1 identified by Singh et al., Mol Cell Biol 2006; doi:10.1128/MCB.26.4.1333-1346.2006). ASO efficacy depends on target accessibility:
- **FIRM stems** (high P) indicate stable secondary structure that may reduce binding.
- **Unpaired regions or FLOPPY stems** (low P) suggest accessible sites.

FoldTrust tiers could triage ASO candidate sites by accessibility. However:
- *Limitation:* ISS-N1 structure in vivo may differ from naked RNA predictions due to protein binding (hnRNP A1, other RBPs).
- *Validation needed:* Correlate tier assignments with ASO binding affinity or knockdown efficacy in cell assays.

**Citation:** Hua Y, Vickers TA, Okunola HL, Bennett CF, Krainer AR. Antisense Masking of an hnRNP A1/A2 Intronic Splicing Silencer Corrects SMN2 Splicing in Transgenic Mice. Am J Hum Genet 2008;82:834-848. doi:10.1016/j.ajhg.2008.01.014

---

## (b) SARS-CoV-2 FSE: Frameshift-Stimulating Element as Small-Molecule Target

**Disease:** COVID-19  
**Target:** Frameshift-stimulating element (NC_045512.2:13462-13542)  
**Therapeutic strategy:** Small molecules disrupting −1 ribosomal frameshifting

**FoldTrust result:** 5 MFE stems: 2 FIRM (mean p 0.946, 0.981) and 3 SOFT (0.730, 0.847, 0.564); ViennaRNA cannot represent the FSE pseudoknot.

**Drug-discovery question:** Can small molecules bind FIRM stems to alter frameshifting efficiency?

The FSE directs −1 programmed ribosomal frameshifting (PRF) required for ORF1ab translation. Disrupting FSE structure reduces PRF efficiency and viral replication:
- Merafloxacin (a fluoroquinolone) was reported to inhibit SARS-CoV-2 −1 PRF and restrict replication in cell culture (Sun et al., PNAS 2021; doi:10.1073/pnas.2023051118). MTDB reduced SARS-CoV-2 −1 PRF in reporter assays (Kelly et al., J Biol Chem 2020; doi:10.1074/jbc.AC120.013449; Neupane et al., J Mol Biol 2020; doi:10.1016/j.jmb.2020.09.006).

**Layer 4 findings:**
- SHAPE reactivity (unpaired probability) correlates with predictions (ρ=0.29–0.55).
- MEA F1 0.83–0.98 vs unconstrained; MFE F1 0.40–0.96 (Pyle bp distance 30) (benchmarks/outputs/layer4_shape/shape_directed_folding.csv).
- FSE correlation is not exceptional genome-wide (43rd–82nd percentile; empirical p 0.18–0.57; layer4_shape_results.json).

FoldTrust FIRM tiers help identify stable stem-loop targets for fragment-based drug design or high-throughput screening. However:
- *Limitation:* Thermodynamic model only; does not predict ligand binding pockets.
- *Validation needed:* Crystallographic or cryo-EM structure of FSE (available: PDB 6XRZ, Zhang et al., Nat Struct Mol Biol 2021; doi:10.1038/s41594-021-00653-y) for structure-based design.

**Citations:**
- Sun Y, Abriola L, Niederer RO, Pedersen SF, Alfajaro MM, Silva Monteiro V, Wilen CB, Ho YC, Gilbert WV, Surovtseva YV, Lindenbach BD, Guo JU. Restriction of SARS-CoV-2 replication by targeting programmed −1 ribosomal frameshifting. Proc Natl Acad Sci USA 2021;118:e2023051118. doi:10.1073/pnas.2023051118
- Kelly JA, Olson AN, Neupane K, Munshi S, San Emeterio J, Pollack L, Woodside MT, Dinman JD. Structural and functional conservation of the programmed −1 ribosomal frameshift signal of SARS coronavirus 2 (SARS-CoV-2). J Biol Chem 2020;295:10741-10748. doi:10.1074/jbc.AC120.013449
- Neupane K, Munshi S, Zhao M, Ritchie DB, Ileperuma SM, Woodside MT. Anti-Frameshifting Ligand Active against SARS Coronavirus-2 Is Resistant to Natural Mutations of the Frameshift-Stimulatory Pseudoknot. J Mol Biol 2020;432:5843-5847. doi:10.1016/j.jmb.2020.09.006

---

## (c) HCV IRES Domain II: Small-Molecule/Ligand Binding Site

**Disease:** Hepatitis C virus infection  
**Target:** Internal ribosome entry site (IRES) domain II (AF009606.1:44-118)  
**Therapeutic strategy:** Small molecules or peptides disrupting IRES-mediated translation

**FoldTrust result:** 2 FIRM / 3 SOFT stems, verdict NEED PROBING (case_metrics.csv); every FIRM and SOFT stem is lost when 25–100 nt of genomic flank is added (Layer 5 window_context_stem_retention.csv).

**Drug-discovery question:** Can ligands bind domain II to inhibit IRES function?

The HCV IRES mediates cap-independent translation initiation, essential for viral protein synthesis. Domain II contributes to ribosome recruitment:
- Benzimidazole ligands bind the domain IIa internal loop and inhibit IRES function (Seth et al., J Med Chem 2005; doi:10.1021/jm050815o; Parsons et al., Nat Chem Biol 2009; doi:10.1038/nchembio.217). Note the ligand site is a loop, not a FIRM stem.
- Domain II adopts a conserved structure (PDB 1P5P; Lukavsky et al., Nat Struct Biol 2003; doi:10.1038/nsb1004).

FoldTrust can prioritize IRES domains for ligand screening. However:
- *Limitation:* Tertiary structure (3D folds, long-range interactions) not predicted by secondary-structure models.
- *Validation needed:* NMR or crystallography to identify binding pockets.

**Citations:**
- Kieft JS, Zhou K, Jubin R, Doudna JA. Mechanism of ribosome recruitment by hepatitis C IRES RNA. RNA 2001;7:194-206. doi:10.1017/s1355838201001790
- Seth PP, Miyaji A, Jefferson EA, Sannes-Lowery KA, Osgood SA, Propp SS, Ranken R, Massire C, Sampath R, Ecker DJ, Swayze EE, Griffey RH. SAR by MS: Discovery of a New Class of RNA-Binding Small Molecules for the Hepatitis C Virus: Internal Ribosome Entry Site IIA Subdomain. J Med Chem 2005;48:7099-7102. doi:10.1021/jm050815o
- Parsons J, Castaldi MP, Dutta S, Dibrov SM, Wyles DL, Hermann T. Conformational inhibition of the hepatitis C virus internal ribosome entry site RNA. Nat Chem Biol 2009;5:823-825. doi:10.1038/nchembio.217
- Lukavsky PJ, Kim I, Otto GA, Puglisi JD. Structure of HCV IRES domain II determined by NMR. Nat Struct Biol 2003;10:1033-1038. doi:10.1038/nsb1004

---

## (d) CFTR 5'UTR: Translation Regulatory Structure

**Disease:** Cystic fibrosis  
**Target:** 5' untranslated region (NM_000492.4:1-200)  
**Therapeutic strategy:** Modulate translation efficiency or mRNA stability

**FoldTrust result:** 3 FIRM / 1 SOFT / 6 FLOPPY stems, verdict REDESIGN (case_metrics.csv); the 3 FIRM stems are retained under every temperature, parameter set and flank (Layer 5). Window = 5′UTR nt 1–70 + first 130 nt CDS.

**Drug-discovery question:** Do 5'UTR structures regulate CFTR translation, and can they be targeted?

CFTR 5'UTR structure influences translation initiation. 5′UTR variants can alter CFTR translation (Lukowski et al., Hum Mutat 2011; doi:10.1002/humu.21545):
- **FIRM stems** in the 5'UTR may create barriers to ribosome scanning.
- Small molecules or ASOs that destabilize inhibitory structures could enhance CFTR translation.

However:
- *Limitation:* 5'UTR effects depend on cellular context (ribosomes, translation factors, RNA-binding proteins).
- *Clinical relevance uncertain:* Current CF therapies (ivacaftor, lumacaftor) target CFTR protein directly; 5'UTR-based strategies remain experimental.

**Citation:** Lukowski SW, Bombieri C, Trezise AEO. Disrupted posttranscriptional regulation of the cystic fibrosis transmembrane conductance regulator (CFTR) by a 5′UTR mutation is associated with a CFTR-related disease. Hum Mutat 2011;32:E2266-E2282. doi:10.1002/humu.21545

---

## (e) MAPT Exon 10: Splicing Regulatory Hairpin

**Disease:** Frontotemporal dementia (FTD) and other tauopathies  
**Target:** MAPT exon 10 regulatory stem-loop (NG_007398.2:120818-121000)  
**Therapeutic strategy:** Antisense oligonucleotides or small molecules modulating splicing

**FoldTrust result:** 5 FIRM / 2 SOFT / 4 FLOPPY stems, verdict REDESIGN (case_metrics.csv); FIRM retention 0.29 at 25-nt flanks, 1.00 at 50/100 nt (Layer 5).

**Drug-discovery question:** Can disrupting or stabilizing the hairpin alter exon 10 inclusion?

MAPT exon 10 encodes a microtubule-binding repeat. Its inclusion is regulated by a stem-loop structure at the 5' splice site:
- FTDP-17 mutations in the hairpin (e.g., S305S, +3, +14, +16) destabilize it and increase exon 10 inclusion (excess 4R tau) (Varani et al., PNAS 1999; doi:10.1073/pnas.96.14.8229; Hutton et al., Nature 1998; doi:10.1038/31508). N279K acts through an exonic splicing enhancer and also increases inclusion (D'Souza et al., PNAS 1999; doi:10.1073/pnas.96.10.5598).
- ASOs shifting exon 10 splicing altered 4R tau in mice (Schoch et al., Neuron 2016; doi:10.1016/j.neuron.2016.04.042); small molecules that bind the splicing-regulatory hairpin have been reported (Chen et al., JACS 2020; doi:10.1021/jacs.0c00768).

FoldTrust FIRM tiers identify the stable hairpin as a high-confidence target. However:
- *Limitation:* Splicing also depends on SR proteins, hnRNPs, and cis elements beyond the hairpin.
- *Validation needed:* Cell-based splicing assays to correlate hairpin stability with exon 10 inclusion.

**Citations:**
- Varani L, Hasegawa M, Spillantini MG, Smith MJ, Murrell JR, Ghetti B, Klug A, Goedert M, Varani G. Structure of tau exon 10 splicing regulatory element RNA and destabilization by mutations of frontotemporal dementia and parkinsonism linked to chromosome 17. Proc Natl Acad Sci USA 1999;96:8229-8234. doi:10.1073/pnas.96.14.8229
- Hutton M, Lendon CL, Rizzu P, Baker M, Froelich S, Houlden H, Pickering-Brown S, Chakraverty S, Isaacs A, Grover A, Hackett J, Adamson J, Lincoln S, Dickson D, Davies P, Petersen RC, Stevens M, de Graaff E, Wauters E, van Baren J, Hillebrand M, Joosse M, Kwon JM, Nowotny P, Che LK, Norton J, Morris JC, Reed LA, Trojanowski J, Basun H, Lannfelt L, Neystat M, Fahn S, Dark F, Tannenberg T, Dodd PR, Hayward N, Kwok JBJ, Schofield PR, Andreadis A, Snowden J, Craufurd D, Neary D, Owen F, Oostra BA, Hardy J, Goate A, van Swieten J, Mann D, Lynch T, Heutink P. Association of missense and 5'-splice-site mutations in tau with the inherited dementia FTDP-17. Nature 1998;393:702-705. doi:10.1038/31508
- D'Souza I, Poorkaj P, Hong M, Nochlin D, Lee VM, Bird TD, Schellenberg GD. Missense and silent tau gene mutations cause frontotemporal dementia with parkinsonism-chromosome 17 type, by affecting multiple alternative RNA splicing regulatory elements. Proc Natl Acad Sci USA 1999;96:5598-5603. doi:10.1073/pnas.96.10.5598
- Schoch KM, DeVos SL, Miller RL, Chun SJ, Norrbom M, Wozniak DF, Dawson HN, Bennett CF, Rigo F, Miller TM. Increased 4R-Tau Induces Pathological Changes in a Human-Tau Mouse Model. Neuron 2016;90:941-947. doi:10.1016/j.neuron.2016.04.042
- Chen JL, Zhang P, Abe M, Talukdar I, Wang F, Li S, Varani G, Hecht SM, Angell L, Wan Y, Childs-Disney JL, Disney MD. Design, Optimization, and Study of Small Molecules That Target Tau Pre-mRNA and Affect Splicing. J Am Chem Soc 2020;142:8706-8727. doi:10.1021/jacs.0c00768

---

## Caveats and Limitations

1. **Hypothesis-level only:** These are potential applications, not validated workflows. No cell-based or animal model data support FoldTrust tier use in drug discovery.

2. **In vitro vs. in vivo:** FoldTrust predicts naked RNA structure. Cellular context (proteins, ions, crowding) alters folding.

3. **No binding prediction:** Tiers indicate structural confidence, not ligand binding sites or affinities. Requires 3D structure (NMR, cryo-EM, X-ray) for rational design.

4. **FLOPPY tier uninformative:** Low-probability stems should not guide targeting decisions.

5. **Validation gap:** No published studies correlate FoldTrust tiers with ASO efficacy, ligand binding, or functional outcomes.

---

## Conclusion

FIRM/SOFT/FLOPPY tiers are a triage signal for which predicted structure is ensemble-supported: FIRM regions are candidates for structure-dependent ligand hypotheses, while less-structured regions are candidates for ASO accessibility hypotheses. Neither is a binding or efficacy prediction.

However, experimental validation is essential. Tier predictions should complement, not replace, cell-based assays, structural biology, and clinical data.

---

**DOIs checked against Crossref (title + first author); see benchmarks/outputs/doi_check.txt.**
