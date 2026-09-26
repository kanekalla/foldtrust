#!/usr/bin/env python3
"""Check ViennaRNA Python API parameter loaders and whether they change MFE energies.
Test sequences are taken from bundle files (no hand-typed sequences):
  1) first record of rfam/rfam_seed_ss.jsonl (RF00005 tRNA)
  2) SARS-CoV-2 FSE region NC_045512.2:13462-13542 (1-based inclusive) from genome/NC_045512.2.fasta
Each measurement runs in a fresh subprocess AND builds the fold_compound with a freshly created
RNA.md() after loading parameters. PITFALL found: in ViennaRNA 2.7.2, RNA.fold(seq) and
RNA.fold_compound(seq) without an explicit md can keep using the previously cached parameter set
after params_load*(), so naive in-process comparisons report identical energies.
"""

import json
import os
import subprocess
import sys

here = os.path.dirname(os.path.abspath(__file__))
root = os.path.dirname(here)


def one(loader, seq):
    code = (
        "import RNA,sys,json\n"
        f"L={loader!r}\n"
        "rc = RNA.params_load(L) if L.endswith('.par') else getattr(RNA, L)()\n"
        "md = RNA.md(); fc = RNA.fold_compound(sys.argv[1], md); st, e = fc.mfe()\n"
        "print(json.dumps({'rc': rc, 'structure': st, 'mfe_kcal_mol': round(e, 2)}))\n"
    )
    r = subprocess.run(
        [sys.executable, "-c", code, seq], capture_output=True, text=True, check=True
    )
    return json.loads(r.stdout)


def main():
    import RNA

    rec = json.loads(open(os.path.join(root, "rfam/rfam_seed_ss.jsonl")).readline())
    g = "".join(
        line.strip()
        for line in open(os.path.join(root, "genome/NC_045512.2.fasta"))
        if not line.startswith(">")
    )
    tests = {
        f"{rec['family']}:{rec['seq_id']}": rec["sequence"],
        "NC_045512.2:13462-13542 (FSE)": g[13461:13542].replace("T", "U"),
    }
    out = {
        "ViennaRNA_version": RNA.__version__,
        "api_has": {
            n: hasattr(RNA, n)
            for n in [
                "params_load_RNA_Turner2004",
                "params_load_RNA_Andronescu2007",
                "params_load_RNA_Langdon2018",
                "params_load",
            ]
        },
        "method": "fresh process per measurement; params loaded, then RNA.md() created, then fold_compound(seq, md).mfe()",
        "results": {},
    }
    loaders = [
        "params_load_RNA_Turner2004",
        "params_load_RNA_Andronescu2007",
        "params_load_RNA_Langdon2018",
    ] + [
        os.path.join(root, "vienna_params", f)
        for f in ["rna_turner2004.par", "rna_andronescu2007.par", "rna_langdon2018.par"]
    ]
    for name, s in tests.items():
        r = {"length": len(s), "sequence": s}
        for L in loaders:
            r[os.path.basename(L)] = one(L, s)
        out["results"][name] = r
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
