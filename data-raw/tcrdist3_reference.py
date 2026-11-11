"""Compute tcrdist3 reference distances (query x clone, alpha and beta) for the
tcrhex regression test.

Usage: python data-raw/tcrdist3_reference.py queries.tsv clones.tsv out.tsv
Both inputs: tab-separated with columns id va ca vb cb (V genes without allele;
'*01' is appended as tcrdist3 expects). Output: qid cid d_a d_b for ALL pairs.
"""
import sys
import pandas as pd
from tcrdist.repertoire import TCRrep

def rep(df):
    cell = pd.DataFrame({
        "count": 1,
        "v_a_gene": df["va"] + "*01", "cdr3_a_aa": df["ca"],
        "v_b_gene": df["vb"] + "*01", "cdr3_b_aa": df["cb"],
        "id": df["id"].astype(str)})
    return TCRrep(cell_df=cell, organism="human", chains=["alpha", "beta"],
                  deduplicate=False, compute_distances=False,
                  db_file="alphabeta_gammadelta_db.tsv")

q = pd.read_csv(sys.argv[1], sep="\t", dtype=str)
c = pd.read_csv(sys.argv[2], sep="\t", dtype=str)
tq = rep(q)
tc = rep(c)
tc.compute_rect_distances(df=tq.clone_df, df2=tc.clone_df)
qa = tq.clone_df["id"].tolist()
ca = tc.clone_df["id"].tolist()
ra, rb = tc.rw_alpha, tc.rw_beta
with open(sys.argv[3], "w") as fh:
    fh.write("qid\tcid\td_a\td_b\n")
    for i, qi in enumerate(qa):
        for j, cj in enumerate(ca):
            fh.write(f"{qi}\t{cj}\t{int(ra[i, j])}\t{int(rb[i, j])}\n")
print(len(qa), "x", len(ca))
