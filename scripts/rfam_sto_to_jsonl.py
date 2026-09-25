#!/usr/bin/env python3
"""Convert Rfam seed Stockholm alignments (with #=GC SS_cons, WUSS notation) into a JSONL
of individual seed sequences with the consensus structure projected onto each sequence.

Rules (all deterministic, no manual edits):
  * Stockholm blocks are concatenated per sequence name (multi-block safe).
  * WUSS -> pairs:
      - mode "wuss" (default, WUSS-correct): (), <>, [], {} are all *nested* helix brackets
        in WUSS and are parsed as base pairs with one shared stack (type-matched).
      - mode "round_angle_only" (literal reading of the task text): only () and <> are
        parsed as pairs; [] and {} columns are treated as unpaired.
      - In both modes pseudoknot letters (Aa, Bb, ...) are stripped (treated unpaired), and all
        other WUSS symbols (, . : _ - ~ etc.) are unpaired.
  * Projection onto each sequence: gap columns (. - _ ~) are removed; a consensus pair is kept
    only if both columns are residues in that sequence and the residue pair is canonical
    Watson-Crick or GU wobble (AU UA GC CG GU UG). T is converted to U, sequence is upper-cased.
  * Sequences containing non-ACGU residues after conversion (e.g. N, R, Y) are skipped.
  * Sequences with length > 400 nt are skipped.
  * Per family: sort remaining records by seq_id (plain string sort) and keep the first 25.
  * Every record is verified: len(sequence) == len(structure), only ACGU, only '().',
    brackets balanced, every pair canonical/wobble.
Usage: rfam_sto_to_jsonl.py OUT.jsonl [--mode wuss|round_angle_only] STO_FILE...
"""
import json, os, sys, re
from collections import OrderedDict

MAX_LEN = 400
MAX_PER_FAMILY = 25
CANON = {("A","U"),("U","A"),("G","C"),("C","G"),("G","U"),("U","G")}
GAP = set(".-_~")
URL = "https://rfam.org/family/{acc}/alignment/stockholm"

def parse_stockholm(path):
    seqs = OrderedDict(); ss = []; gf = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip() or line.startswith("# STOCKHOLM") or line.startswith("//"):
                continue
            if line.startswith("#=GC SS_cons"):
                ss.append(line.split(None, 2)[2].strip())
            elif line.startswith("#=GF "):
                p = line.split(None, 2)
                if len(p) == 3 and p[1] in ("AC", "ID", "DE"):
                    gf.setdefault(p[1], p[2].strip())
            elif line.startswith("#"):
                continue
            else:
                name, s = line.split(None, 1)
                seqs[name] = seqs.get(name, "") + s.strip()
    ss = "".join(ss)
    for n, s in seqs.items():
        if len(s) != len(ss):
            raise ValueError(f"{path}: {n} aligned length {len(s)} != SS_cons length {len(ss)}")
    return gf, seqs, ss

def wuss_pairs(ss, mode):
    opens = {"(": ")", "<": ">", "[": "]", "{": "}"} if mode == "wuss" else {"(": ")", "<": ">"}
    closes = {v: k for k, v in opens.items()}
    stack = []; pairs = {}
    for i, c in enumerate(ss):
        if c in opens:
            stack.append((c, i))
        elif c in closes:
            if not stack or stack[-1][0] != closes[c]:
                raise ValueError(f"unbalanced/crossing WUSS bracket at column {i}")
            _, j = stack.pop(); pairs[j] = i; pairs[i] = j
        # letters (pseudoknots) and all other symbols -> unpaired
    if stack:
        raise ValueError("unbalanced WUSS: unclosed brackets")
    return pairs

def project(aln, pairs):
    cols = [i for i, c in enumerate(aln) if c not in GAP]
    col2pos = {c: k for k, c in enumerate(cols)}
    seq = "".join(aln[c] for c in cols).upper().replace("T", "U")
    db = ["."] * len(seq)
    for i, j in pairs.items():
        if i < j and i in col2pos and j in col2pos:
            a, b = col2pos[i], col2pos[j]
            if (seq[a], seq[b]) in CANON:
                db[a] = "("; db[b] = ")"
    return seq, "".join(db)

def verify(rec):
    s, d = rec["sequence"], rec["structure"]
    assert len(s) == len(d), rec["seq_id"]
    assert set(s) <= set("ACGU"), rec["seq_id"]
    assert set(d) <= set("()."), rec["seq_id"]
    st = []
    for k, c in enumerate(d):
        if c == "(":
            st.append(k)
        elif c == ")":
            assert st, f"unbalanced {rec['seq_id']}"
            j = st.pop(); assert (s[j], s[k]) in CANON, rec["seq_id"]
    assert not st, f"unbalanced {rec['seq_id']}"

def main():
    args = sys.argv[1:]
    out = args.pop(0); mode = "wuss"
    if args and args[0] == "--mode":
        args.pop(0); mode = args.pop(0)
    assert mode in ("wuss", "round_angle_only")
    stats = OrderedDict(); n = 0
    with open(out, "w") as fo:
        for path in args:
            gf, seqs, ss = parse_stockholm(path)
            acc = gf.get("AC") or os.path.basename(path).split(".")[0]
            pairs = wuss_pairs(ss, mode)
            recs = []; skipped = {"non_ACGU": 0, "too_long": 0}
            for name in sorted(seqs):
                seq, db = project(seqs[name], pairs)
                if not set(seq) <= set("ACGU"):
                    skipped["non_ACGU"] += 1; continue
                if len(seq) > MAX_LEN:
                    skipped["too_long"] += 1; continue
                recs.append({"family": acc, "seq_id": name, "sequence": seq, "structure": db,
                             "source_url": URL.format(acc=acc)})
            recs = recs[:MAX_PER_FAMILY]
            for r in recs:
                verify(r); fo.write(json.dumps(r) + "\n"); n += 1
            stats[acc] = {"id": gf.get("ID"), "seed_seqs": len(seqs), "kept": len(recs), **skipped,
                          "mean_pairs": round(sum(r["structure"].count("(") for r in recs) / max(1, len(recs)), 1)}
    json.dump({"mode": mode, "total_records": n, "families": stats}, sys.stderr, indent=1)
    sys.stderr.write("\n")

if __name__ == "__main__":
    main()
