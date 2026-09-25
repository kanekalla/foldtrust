#!/usr/bin/env python3
"""Convert ArchiveII .ct files (read directly from the original archiveII.tar.gz) to JSONL.
Fields: family (filename prefix before first '_'), name (file stem), sequence (upper-case, T->U,
as in the ct file; non-ACGU letters kept), pairs (1-based [i,j], i<j, exactly as in the ct file,
including non-canonical pairs), structure (dot-bracket; pseudoknotted pairs are assigned to
successive bracket levels () [] {} <> greedily in order of the 5' partner), n_pseudoknot_pairs,
source_url, source_member.  Only the first structure of each ct file is used (count recorded).
Verification per record: pair table symmetric, partners in range, no self pairs, structure length ==
sequence length, brackets balanced per level.
Usage: archiveII_ct_to_jsonl.py archiveII.tar.gz OUT.jsonl SOURCE_URL
"""
import json, sys, tarfile, os, collections

LEVELS = ["()", "[]", "{}", "<>"]

def parse_ct(text):
    lines = [l for l in text.splitlines() if l.strip()]
    n = int(lines[0].split()[0]); rows = lines[1:1 + n]
    n_structs = sum(1 for l in lines if len(l.split()) >= 1 and not l.split()[0].isdigit()) # unused
    seq = []; partner = [0] * (n + 1)
    for k, l in enumerate(rows, 1):
        f = l.split(); assert int(f[0]) == k, f"bad index line {k}"
        seq.append(f[1]); partner[k] = int(f[4])
    extra = (len(lines) - 1) // (n + 1) if n else 0
    return "".join(seq).upper().replace("T", "U"), partner, n, len(lines) > n + 1

def to_dotbracket(n, pairs):
    db = ["."] * n; levels = [[] for _ in LEVELS]; unplaced = 0
    for i, j in sorted(pairs):
        for L, lv in enumerate(levels):
            if all(not (a < i < b < j or i < a < j < b) for a, b in lv):
                lv.append((i, j)); db[i-1], db[j-1] = LEVELS[L][0], LEVELS[L][1]; break
        else:
            unplaced += 1
    return "".join(db), sum(len(l) for l in levels[1:]), unplaced

def main(tgz, out, url):
    stats = collections.Counter(); multi = 0
    with tarfile.open(tgz) as tf, open(out, "w") as fo:
        members = sorted((m for m in tf.getmembers() if m.isfile() and m.name.endswith(".ct")), key=lambda m: m.name)
        for m in members:
            text = tf.extractfile(m).read().decode("latin-1")
            seq, partner, n, has_more = parse_ct(text)
            multi += has_more
            assert len(seq) == n
            pairs = []
            for i in range(1, n + 1):
                j = partner[i]
                if j:
                    assert 1 <= j <= n and j != i and partner[j] == i, f"{m.name}: asymmetric pair {i}-{j}"
                    if i < j: pairs.append([i, j])
            db, npk, unplaced = to_dotbracket(n, pairs)
            assert len(db) == n and unplaced == 0, m.name
            for o, c in LEVELS:
                d = 0
                for ch in db:
                    d += (ch == o) - (ch == c); assert d >= 0
                assert d == 0
            stem = os.path.basename(m.name)[:-3]
            fam = stem.split("_")[0]
            fo.write(json.dumps({"family": fam, "name": stem, "length": n, "sequence": seq, "pairs": pairs,
                                 "structure": db, "n_pseudoknot_pairs": npk, "source_url": url,
                                 "source_member": m.name}) + "\n")
            stats[fam] += 1
    print(json.dumps({"total": sum(stats.values()), "per_family": dict(sorted(stats.items())),
                      "ct_files_with_more_than_one_structure": multi}, indent=1))

if __name__ == "__main__":
    main(*sys.argv[1:4])
