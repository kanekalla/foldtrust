#!/usr/bin/env python3
"""Locate the SARS-CoV-2 ORF1ab -1 frameshift element in NC_045512.2 (files in genome/).
Slippery site found by exact search for UUUAAAC (as TTTAAAC in DNA FASTA); stem-loops and the
ribosomal_slippage CDS join are read from the NCBI GenBank flat file features. 1-based inclusive."""
import json, re, os
root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
g = "".join(l.strip() for l in open(os.path.join(root, "genome/NC_045512.2.fasta")) if not l.startswith(">"))
hits = [m.start() + 1 for m in re.finditer("(?=TTTAAAC)", g)]
near = [h for h in hits if 13000 <= h <= 14000]
assert len(near) == 1
s = near[0]; e = s + 6
gb = open(os.path.join(root, "genome/NC_045512.2.gb")).read()
stems = re.findall(r"stem_loop\s+(\d+)\.\.(\d+)\s+/gene=\"ORF1ab\".*?/function=\"([^\"]+)\"", gb, re.S)
stems = [(int(a), int(b), " ".join(f.split())) for a, b, f in stems if 13000 < int(a) < 14000]
cds = re.search(r"CDS\s+(join\(266\.\.13468,13468\.\.21555\))", gb)
fse_end = max(b for _, b, _ in stems)
out = {"reference": "NC_045512.2", "genome_length": len(g), "coordinate_system": "1-based, inclusive",
       "slippery_site": {"start": s, "end": e, "sequence_RNA": g[s-1:e].replace("T", "U")},
       "all_UUUAAAC_occurrences_in_genome": hits,
       "genbank_stem_loops": [{"start": a, "end": b, "function": f, "sequence_RNA": g[a-1:b].replace("T", "U")} for a, b, f in stems],
       "genbank_ORF1ab_CDS_ribosomal_slippage": cds.group(1) if cds else None,
       "fse_span_slippery_site_to_end_of_stem_loop_2": {"start": s, "end": fse_end,
            "length": fse_end - s + 1, "sequence_RNA": g[s-1:fse_end].replace("T", "U")}}
print(json.dumps(out, indent=1))
