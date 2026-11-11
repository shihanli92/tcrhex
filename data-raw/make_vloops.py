"""Build inst/extdata/tcrdist_vloops_human.tsv from tcrdist3's alphabeta_gammadelta_db.tsv.

tcrdist3 (MIT License, Copyright (c) 2020 Koshlan Mayer-Blackwell) ships the
IMGT-aligned CDR1 / CDR2 / CDR2.5 (pMHC-facing) loops of every V allele. tcrdist3
looks loops up by allele (e.g. "TRAV27*01"); we key by gene and use the *01 allele
(or the lowest-numbered allele with loops if *01 is absent), recording which.

Includes the gamma/delta V genes the db files under chains A/B: tcrdist3 looks
loops up by gene name whatever the chain (e.g. TRDV1 rearranged into TRA).

Usage: python data-raw/make_vloops.py <path/to/alphabeta_gammadelta_db.tsv> inst/extdata/tcrdist_vloops_<organism>.tsv [human|mouse]
"""
import csv, sys, re
src, out = sys.argv[1], sys.argv[2]
organism = sys.argv[3] if len(sys.argv) > 3 else "human"
rows = {}
with open(src) as fh:
    for l in csv.DictReader(fh, delimiter="\t"):
        if l["organism"] != organism or l["region"] != "V" or l["chain"] not in ("A", "B"):
            continue
        if not l["cdrs"]:
            continue
        cdrs = l["cdrs"].split(";")
        gene, allele = l["id"].split("*")
        rows.setdefault(gene, []).append((int(allele), l["id"], l["chain"], cdrs[0], cdrs[1], cdrs[2]))
with open(out, "w") as fh:
    fh.write("gene\tchain\tallele\tcdr1\tcdr2\tcdr25\tother_alleles_differ\n")
    for gene in sorted(rows):
        alle = sorted(rows[gene])
        a = alle[0]
        differ = any(x[3:6] != a[3:6] for x in alle[1:])
        fh.write("\t".join([gene, "TR" + a[2], a[1], a[3], a[4], a[5], str(differ).upper()]) + "\n")
print(len(rows), "genes")
