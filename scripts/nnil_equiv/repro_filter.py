#!/usr/bin/env python
"""
Reproduce Holland's File_S10 filter on File S1 exactly, and diff the result against
(a) the published File S18 marker/line set and (b) our current geno_recoded, to
locate the 35-state residual. Faithful transcription of File_S10:
 merge File S9 (inner) -> clean names -> drop duplicate names (keep first) ->
 recode (x-1)*-2 -> keep markers 0<maf<=0.05 -> drop markers het>=0.02 ->
 drop lines het>=0.02.
  ~/anaconda3/envs/nilhmm/bin/python scripts/nnil_equiv/repro_filter.py
"""
import os
import numpy as np
import pandas as pd

ROOT = os.path.expanduser("~/repos/zealhmm")
ZH   = os.path.join(ROOT, "data/zhong2025")
NNIL = os.path.join(ROOT, "agent/nNIL")
OUT  = os.path.join(ROOT, "data/nnil_equiv")


def log(m): print(f"[repro_filter] {m}", flush=True)


# --- File_S10, faithfully ---------------------------------------------------
geno = pd.read_table(os.path.join(ZH, "File_S01.nNIL_raw_SNPs_bgi_id_miss20.txt"),
                     sep="\t", skiprows=1, header="infer")
marker_names = geno.columns.to_series()[1:]                 # File S1 column order
trans = pd.read_table(os.path.join(NNIL, "File_S09.bgi_nil_id.txt"),
                      header=0, encoding="windows-1252")
genob = pd.merge(trans, geno, how="inner", left_on="bgi_id", right_on="<Marker>")
genob = genob.drop(labels=["bgi_id", "<Marker>"], axis=1)
sample_names = genob.iloc[:, 0]
sample_names = sample_names.str.replace("\xa0", " ").str.replace("Â", "")  # \xa0, Â
dup = sample_names.duplicated()
log(f"merged {genob.shape[0]} samples; duplicate names dropped: {int(dup.sum())}")
genob = genob.loc[~dup, :]; sample_names = sample_names[~dup]
genob = genob.iloc[:, 1:]

genomat = np.multiply(np.subtract(genob.to_numpy(), 1), -2)  # {0,1,2}, NaN
maf = np.multiply(np.nanmean(genomat, axis=0), 0.5)
keepm = ~((maf == 0) | (maf > 0.05))
mnf = marker_names[np.asarray(keepm)]; g = genomat[:, np.asarray(keepm)]
hetm = np.nanmean(g == 1, axis=0)
g = g[:, hetm < 0.02]; mnf = mnf[hetm < 0.02]
hetl = np.nanmean(g == 1, axis=1)
g = g[hetl < 0.02, :]; snf = sample_names[hetl < 0.02]
log(f"File_S10 filtered set: {g.shape[0]} lines x {g.shape[1]} markers")

# --- diff vs File S18 (published) set ---------------------------------------
s18 = pd.read_csv(os.path.join(ZH, "File_S18.nNIL_gbs_HMM_introgressionCalls_full_set.csv"),
                  index_col=0, nrows=0)
s18_markers = set(s18.columns); s18_lines = None
s18idx = pd.read_csv(os.path.join(ZH, "File_S18.nNIL_gbs_HMM_introgressionCalls_full_set.csv"),
                     index_col=0, usecols=[0]).index.astype(str)
mnf_set = set(mnf.tolist()); snf_set = set(snf.astype(str).tolist())
log(f"markers: File_S10={len(mnf_set)}  File S18={len(s18_markers)}  "
    f"S10-only={len(mnf_set - s18_markers)}  S18-only={len(s18_markers - mnf_set)}")
log(f"lines:   File_S10={len(snf_set)}  File S18={len(s18idx)}  "
    f"S10-only={len(snf_set - set(s18idx))}  S18-only={len(set(s18idx) - snf_set)}")

# --- diff genotypes vs our current geno_recoded on the 12 flagged lines ------
gr = pd.read_csv(os.path.join(OUT, "geno_recoded.csv"), index_col=0)
gr.index = gr.index.astype(str)
filt = pd.DataFrame(np.where(np.isnan(g), 3, g).astype(int),
                    index=snf.astype(str).values, columns=mnf.values)
common_l = [l for l in gr.index if l in set(filt.index)]
common_m = [m for m in gr.columns if m in set(filt.columns)]
a = gr.loc[common_l, common_m].to_numpy()
b = filt.loc[common_l, common_m].to_numpy()
d = a != b
log(f"geno_recoded vs File_S10 genotypes on {len(common_l)} lines x {len(common_m)} markers: "
    f"{int(d.sum())} cell differences")
if d.sum():
    rows, cols = np.where(d)
    from collections import Counter
    log(f"  differing lines: {len(set(rows))}; per-line: {sorted(Counter(rows).values(), reverse=True)[:12]}")
    log(f"  my value -> File_S10 value at diffs: {Counter(zip(a[rows,cols], b[rows,cols]))}")
