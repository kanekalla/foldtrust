#!/usr/bin/env python3
"""Generate SHA256SUMS and MANIFEST.md for the bundle (metadata below; hashes/sizes computed)."""

import glob
import hashlib
import os

root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
cache_dir = os.path.join(root, "data", "_cache")
os.chdir(cache_dir)
DAS = "75cb3151e200dd7ea62600db50026ce179181fa3"
VIE = "1ffec79f5e258896160f7362ced8263450f371dc"
C = {
    "archiveII": "Sloma MF, Mathews DH. Exact calculation of loop formation probability identifies folding motifs in RNA secondary structures. RNA 2016;22:1808-1818 (ArchiveII; earlier versions Mathews et al. 1999 JMB 288:911; Bellaousov & Mathews 2010 RNA 16:1870). Original structure sources listed in archiveII/README inside the tarball must also be cited.",
    "mxfold2": "Sato K, Akiyama M, Sakakibara Y. RNA secondary structure prediction using deep learning with thermodynamic integration. Nat Commun 2021;12:941. Dataset DOI 10.5281/zenodo.4430150. Contained sets: TrainSetA/TestSetA/TrainSetB/TestSetB (Rivas E, Lang R, Eddy SR. RNA 2012;18:193-212); bpRNA-1m TR0/VL0/TS0 (Danaee P et al. NAR 2018;46:5381-5394; split from Singh J et al. Nat Commun 2019); RNAStrAlign (Tan Z et al. Bioinformatics 2017;33:3698); archiveII (Sloma & Mathews 2016); bpRNA-new (Sato et al. 2021).",
    "spot": "Singh J, Hanson J, Paliwal K, Zhou Y. RNA secondary structure prediction using an ensemble of two-dimensional deep neural networks and transfer learning. Nat Commun 2019;10:5407.",
    "rfam": "Ontiveros-Palacios N, Cooke E, Nawrocki EP, Triebel S, Marz M, Rivas E, Griffiths-Jones S, Petrov AI, Bateman A, Sweeney B. Rfam 15: RNA families database in 2025. Nucleic Acids Res 2025;53(D1):D258-D267. Rfam release current at download = 15.1 (per FTP CURRENT/README). Licence CC0.",
    "das_repo": f"DasLab/SARS_CoV-2_shape_comparison @ {DAS} (README). Normalization per descriptions.txt.",
    "incarnato": "Manfredonia I, et al. Genome-wide mapping of SARS-CoV-2 RNA structures identifies therapeutically-relevant elements. Nucleic Acids Res 2020;48:12436-12452 (Incarnato lab).",
    "pyle": "Huston NC, et al. Comprehensive in vivo secondary structure of the SARS-CoV-2 genome reveals novel regulatory motifs and mechanisms. Mol Cell 2021;81:584-598 (Pyle lab).",
    "zhang": "Sun L, et al. In vivo structural characterization of the SARS-CoV-2 RNA genome identifies host proteins vulnerable to repurposed drugs. Cell 2021;184:1865-1883 (Zhang lab).",
    "rangan": "Rangan R, et al. RNA genome conservation and secondary structure in SARS-CoV-2 and SARS-related viruses: a first look. RNA 2020;26:937-959.",
    "ncbi": "NCBI RefSeq NC_045512.2 (SARS-CoV-2 Wuhan-Hu-1). Wu F, et al. A new coronavirus associated with a human respiratory disease in China. Nature 2020;579:265-269.",
    "vienna": "Lorenz R, et al. ViennaRNA Package 2.0. Algorithms Mol Biol 2011;6:26.",
    "turner": "Turner 2004 parameters: Mathews DH, et al. PNAS 2004;101:7287-7292 (NNDB, Turner & Mathews NAR 2010;38:D280).",
    "andronescu": "Andronescu M, Condon A, Hoos HH, Mathews DH, Murphy KP. Efficient parameter estimation for RNA secondary structure prediction. Bioinformatics 2007;23:i19-i28.",
    "langdon": "Langdon WB, Petke J, Lorenz R. Evolving better RNAfold structure prediction. EuroGP 2018, LNCS 10781:220-236.",
}
RAW = f"https://raw.githubusercontent.com/DasLab/SARS_CoV-2_shape_comparison/{DAS}/"
M = {}
M["archiveII/archiveII.tar.gz"] = (
    "https://web.archive.org/web/20221017002351id_/https://rna.urmc.rochester.edu/pub/archiveII.tar.gz (Wayback Machine raw capture of the Mathews-lab file; original https://rna.urmc.rochester.edu/pub/archiveII.tar.gz now returns 404). Wayback CDX SHA1 digest MJTVJLRESHZU5ZYN6ZPULE5ASUZPLA2Q verified equal to SHA1 of downloaded bytes.",
    C["archiveII"],
    "Original, unmodified. 3975 .ct files (+ .seq, .lis lists, README).",
)
M["archiveII/archiveII.jsonl.gz"] = (
    "DERIVED from archiveII/archiveII.tar.gz by scripts/archiveII_ct_to_jsonl.py",
    C["archiveII"],
    "3975 records; fields family,name,length,sequence,pairs(1-based, all pairs incl. non-canonical),structure(dot-bracket, PK levels ()[]{}<>),n_pseudoknot_pairs,source_url,source_member. Includes full 16S/23S and their _domainN files.",
)
M["archiveII/archiveII.jsonl.stats.json"] = (
    "DERIVED (converter stdout)",
    C["archiveII"],
    "Per-family record counts.",
)
for k in range(3):
    M[f"tool_data_releases/mxfold2_zenodo_4430150/mxfold2-data.tar.gz.part{k:02d}"] = (
        "https://zenodo.org/api/records/4430150/files/mxfold2-data.tar.gz/content (Zenodo record 4430150)",
        C["mxfold2"],
        "Byte-split (split -b 10000000) of the original 20,360,424-byte mxfold2-data.tar.gz to keep files <20 MB. Reassemble: cat mxfold2-data.tar.gz.part0* > mxfold2-data.tar.gz ; whole-file SHA256 8332469e74c0b2be2140a883f92166488c308e87b662c39eebe771002e1ef51c, MD5 45a4d51f01517dd782216920ccc510b7 (= Zenodo checksum). Contents: archiveII (3966 bpseq), TrainSetA/TestSetA/TrainSetB/TestSetB, bpRNA_dataset-canonicals TR0(10814)/VL0(1300)/TS0(1305), bpRNAnew (5401), RNAStrAlign, .lst files.",
    )
M["tool_data_releases/spot_rna/PDB_dataset.zip"] = (
    "https://www.dropbox.com/s/vnq0k9dg7vynu3q/PDB_dataset.zip?dl=1 (link from SPOT-RNA README)",
    C["spot"],
    "Original. PDB-derived sets TR1/VL1/TS1/TS2 sequences + labels.",
)
fam_names = {
    "RF00005": "tRNA",
    "RF00001": "5S_rRNA",
    "RF00017": "Metazoa_SRP",
    "RF01854": "Bacteria_large_SRP",
    "RF00010": "RNaseP_bact_a",
    "RF00003": "U1",
    "RF00004": "U2",
    "RF00015": "U4",
    "RF00026": "U6",
    "RF00050": "FMN",
    "RF00162": "SAM (SAM-I)",
    "RF00167": "Purine",
    "RF00059": "TPP",
    "RF00174": "Cobalamin",
    "RF00504": "Glycine",
    "RF00380": "Magnesium (ykoK/M-box)",
    "RF00168": "Lysine",
    "RF00234": "glmS",
    "RF00023": "tmRNA",
}
for p in sorted(glob.glob("rfam/raw/*.sto")):
    a = os.path.basename(p)[:-4]
    M[p] = (
        f"https://rfam.org/family/{a}/alignment/stockholm",
        C["rfam"],
        f"Original seed alignment ({fam_names[a]}).",
    )
M["rfam/rfam_seed_ss.jsonl"] = (
    "DERIVED from rfam/raw/*.sto by scripts/rfam_sto_to_jsonl.py (default --mode wuss)",
    C["rfam"],
    "PRIMARY. SS_cons projected per sequence; (), <>, [], {} parsed as nested WUSS pairs; pseudoknot letters stripped; canonical/wobble only; <=400 nt; ACGU-only; first 25 by seq_id. Fields family,seq_id,sequence,structure,source_url.",
)
M["rfam/rfam_seed_ss.round_angle_only.jsonl"] = (
    "DERIVED, same script with --mode round_angle_only",
    C["rfam"],
    "Literal-instruction variant: only () and <> kept as pairs; [] and {} (nested helices in WUSS) treated as unpaired. Same records/order as primary; only structures differ.",
)
M["rfam/rfam_seed_ss.stats.json"] = (
    "DERIVED (converter stderr)",
    C["rfam"],
    "Per-family counts / skip reasons.",
)
M["rfam/rfam_seed_ss.round_angle_only.stats.json"] = (
    "DERIVED (converter stderr)",
    C["rfam"],
    "Per-family counts for variant.",
)
for f, cit, note in [
    (
        "incarnato_invitro_reactivity.csv",
        C["incarnato"],
        "Incarnato lab in vitro SHAPE-MaP; one value per line (29903 lines = 1 per nt of NC_045512.2), no header; nan = no data; range 0-15.16 (not 0-1 normalized)",
    ),
    (
        "incarnato_invivo_reactivity.csv",
        C["incarnato"],
        "Incarnato lab in vivo SHAPE-MaP (same paper; README names only the in vitro set explicitly); 29903 lines, NaN = no data; range 0-8.97",
    ),
    (
        "pyle_reactivity.csv",
        C["pyle"],
        "Pyle lab in vivo SHAPE-MaP; 29903 lines; -999 = no data; range -14.0 to 45.7 (raw, contains negatives)",
    ),
    (
        "zhang_invitro_reactivity.csv",
        C["zhang"],
        "Zhang lab in vitro icSHAPE; 29903 lines; nan = no data; range 0-1",
    ),
    (
        "zhang_invivo_reactivity.csv",
        C["zhang"],
        "Zhang lab in vivo icSHAPE; 29903 lines; nan = no data; range 0-1",
    ),
]:
    M[f"shape/{f}"] = (
        RAW + "SHAPE%20data/" + f,
        cit + " Obtained via " + C["das_repo"],
        note + ". Git blob SHA1 verified against tree at pinned commit.",
    )
M["shape/README.md"] = (RAW + "README.md", C["das_repo"], "Citations for each dataset.")
M["shape/descriptions.txt"] = (
    RAW + "descriptions.txt",
    C["das_repo"],
    "Normalization description.",
)
M["shape/refseq.txt"] = (
    RAW + "refseq.txt",
    C["das_repo"],
    "29903-nt reference; verified identical to NC_045512.2.",
)
M["shape/sarsr_conservation.csv"] = (
    RAW + "sarsr_conservation.csv",
    C["rangan"] + " Via " + C["das_repo"],
    "Conservation track.",
)
for f, cit in [
    ("Incarnato", C["incarnato"]),
    ("pyle", C["pyle"]),
    ("rnaz", C["rangan"]),
    ("zhang", C["zhang"]),
]:
    fn = f"{f}_secondary_structures.csv"
    M[f"shape/secondary_structure_tracks/{fn}"] = (
        RAW + "Secondary%20Structure%20Tracks/" + fn,
        cit + " Via " + C["das_repo"],
        "Proposed structured regions (start,end,seq,dot-bracket,...). Has UTF-8 BOM.",
    )
EU = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_045512.2&rettype={}&retmode=text"
M["genome/NC_045512.2.fasta"] = (EU.format("fasta"), C["ncbi"], "29903 nt.")
M["genome/NC_045512.2.gb"] = (
    EU.format("gb"),
    C["ncbi"],
    "GenBank flat file (features used for FSE stem-loops).",
)
M["genome/fse_coordinates.json"] = (
    "DERIVED by scripts/locate_fse.py",
    C["ncbi"],
    "Slippery site 13462-13468 UUUAAAC; stem-loop 1 13476-13503; stem-loop 2 13488-13542 (GenBank, Rfam RF00507 inference).",
)
for f, cit in [
    ("rna_turner2004.par", C["turner"]),
    ("rna_andronescu2007.par", C["andronescu"]),
    ("rna_langdon2018.par", C["langdon"]),
]:
    M[f"vienna_params/{f}"] = (
        f"https://raw.githubusercontent.com/ViennaRNA/ViennaRNA/{VIE}/misc/{f}",
        cit + " " + C["vienna"],
        "Original (identical to master at download).",
    )
M["vienna_params/param_check_result.json"] = (
    "DERIVED by scripts/vienna_param_check.py with pip ViennaRNA 2.7.2",
    C["vienna"],
    "MFE per parameter set.",
)
M["scripts/manifest_notes.md"] = (
    "written for this bundle",
    "-",
    "Narrative section included in MANIFEST.md.",
)
for p in sorted(glob.glob("scripts/*.py")):
    M[p] = ("written for this bundle", "-", "Script.")


def sha(p):
    h = hashlib.sha256()
    h.update(open(p, "rb").read())
    return h.hexdigest()


allfiles = sorted(os.path.relpath(os.path.join(d, f)) for d, _, fs in os.walk(".") for f in fs)
allfiles = [f for f in allfiles if f not in ("MANIFEST.md", "SHA256SUMS")]
missing = [f for f in allfiles if f not in M]
assert not missing, missing
with open("SHA256SUMS", "w") as fo:
    for f in allfiles:
        fo.write(f"{sha(f)}  {f}\n")
total = sum(os.path.getsize(f) for f in allfiles)
L = []
L.append("# foldtrust-bench-data — public-data bundle for an RNA secondary-structure benchmark\n")
L.append(
    "Built 2026-09-25 (America/Denver) on the shared box at /workspace/foldtrust-bench-data. No git clones were used; every file below was downloaded with curl from the URL given (or derived by a script in scripts/ from such a file). No structure or reactivity value was hand-edited. Verify with `python3 scripts/verify_bundle.py` (checks SHA256SUMS, MXfold2 reassembly, record invariants).\n"
)
L.append(
    f"Total size: {total:,} bytes ({total/1e6:.1f} MB) in {len(allfiles)} files (+MANIFEST.md, SHA256SUMS). Largest file < 20 MB.\n"
)
L.append(open("scripts/manifest_notes.md").read())
L.append(
    "\n## File list\n\n| file | size (bytes) | SHA256 | source URL | citation | notes |\n|---|---:|---|---|---|---|"
)
for f in allfiles:
    u, c, n = M[f]
    L.append(f"| `{f}` | {os.path.getsize(f):,} | `{sha(f)}` | {u} | {c} | {n} |")
open("MANIFEST.md", "w").write("\n".join(L) + "\n")
print("files", len(allfiles), "total bytes", total)
