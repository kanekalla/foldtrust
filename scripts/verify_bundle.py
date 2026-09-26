#!/usr/bin/env python3
"""Verify the bundle: SHA256SUMS, MXfold2 part reassembly, and structure-record invariants."""

import gzip
import hashlib
import json
import os
import sys

root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(root)
ok = True


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


for line in open("SHA256SUMS"):
    d, p = line.split(None, 1)
    p = p.strip()
    if sha(p) != d:
        print("SHA MISMATCH", p)
        ok = False
h = hashlib.sha256()
for k in range(3):
    h.update(
        open(
            f"tool_data_releases/mxfold2_zenodo_4430150/mxfold2-data.tar.gz.part{k:02d}", "rb"
        ).read()
    )
MX = "8332469e74c0b2be2140a883f92166488c308e87b662c39eebe771002e1ef51c"
print("mxfold2 reassembled sha256 ok:", h.hexdigest() == MX)
ok &= h.hexdigest() == MX
CANON = {"AU", "UA", "GC", "CG", "GU", "UG"}


def check_db(s, d, canon_only):
    if len(s) != len(d):
        return "length"
    for o, c in ["()", "[]", "{}", "<>"]:
        st = []
        for i, ch in enumerate(d):
            if ch == o:
                st.append(i)
            elif ch == c:
                if not st:
                    return "unbalanced"
                j = st.pop()
                if canon_only and s[j] + s[i] not in CANON:
                    return "noncanonical"
        if st:
            return "unbalanced"
    return None


for fn in ["rfam/rfam_seed_ss.jsonl", "rfam/rfam_seed_ss.round_angle_only.jsonl"]:
    n = bad = 0
    fam = {}
    for line in open(fn):
        r = json.loads(line)
        n += 1
        fam[r["family"]] = fam.get(r["family"], 0) + 1
        e = check_db(r["sequence"], r["structure"], True)
        if e or set(r["sequence"]) - set("ACGU") or set(r["structure"]) - set("()."):
            bad += 1
    print(fn, "records", n, "bad", bad, "families", len(fam), "max/family", max(fam.values()))
    ok &= bad == 0
n = bad = 0
for line in gzip.open("archiveII/archiveII.jsonl.gz", "rt"):
    r = json.loads(line)
    n += 1
    if check_db(r["sequence"], r["structure"], False):
        bad += 1
print("archiveII records", n, "bad", bad)
ok &= bad == 0
print("ALL OK" if ok else "FAILURES")
sys.exit(0 if ok else 1)
