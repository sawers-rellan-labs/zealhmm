#!/usr/bin/env python
"""
Build the shared nNIL caller input by reproducing Holland's File_S10 filter on the
File S1 raw genotypes -- verified 0-cell-identical to his filtered set (same 888
lines, same 64,025 markers, same values as File S18's set). This supersedes the
earlier "subset to File S18" heuristic (which recovered only 884 lines because it
skipped Holland's name-cleaning); deriving the filter recovers all 888.

File_S10, faithfully: merge File S9 (inner) -> clean sample names (\\xa0, A-hat) ->
drop duplicate names (keep first) -> recode (raw-1)*-2 to {0,1,2} -> keep markers
with 0<maf<=0.05 -> drop markers with het rate >=0.02 -> drop lines with het rate
>=0.02. Missing -> 3 for the caller.

Run with the existing `nilhmm` conda env:
  ~/anaconda3/envs/nilhmm/bin/python scripts/nnil_equiv/01_build_input.py

Outputs (data/nnil_equiv/, gitignored):
  geno_recoded.csv  lines x markers, {0,1,2,3}
  markers.csv       marker, chrom, pos  (File_S10 / File S1 order)
  lines.csv         line names, in row order
  s18_aligned.csv   published File S18 calls, same lines x markers order
  params.json       nir, germ, gert, p, mr, r, f_1, f_2  (File S16 operating point)
"""
import json
import os
import numpy as np
import pandas as pd

ROOT = os.path.expanduser("~/repos/zealhmm")
ZH   = os.path.join(ROOT, "data/zhong2025")
NNIL = os.path.join(ROOT, "agent/nNIL")
OUT  = os.path.join(ROOT, "data/nnil_equiv"); os.makedirs(OUT, exist_ok=True)
S1   = os.path.join(ZH, "File_S01.nNIL_raw_SNPs_bgi_id_miss20.txt")
S18  = os.path.join(ZH, "File_S18.nNIL_gbs_HMM_introgressionCalls_full_set.csv")
S9   = os.path.join(NNIL, "File_S09.bgi_nil_id.txt")
PARAMS = dict(nir=0.9, germ=0.01, gert=0.0001, p=0.9, f_1=0.007813, f_2=0.011179)


def log(m): print(f"[01_build_input] {m}", flush=True)


# ---- File_S10 filter (verbatim logic) --------------------------------------
log("reading File S1 (raw genotypes)...")
geno = pd.read_table(S1, sep="\t", skiprows=1, header="infer")
marker_names = geno.columns.to_series()[1:]                     # File S1 column order
trans = pd.read_table(S9, header=0, encoding="windows-1252")
genob = pd.merge(trans, geno, how="inner", left_on="bgi_id", right_on="<Marker>")
genob = genob.drop(labels=["bgi_id", "<Marker>"], axis=1)
sample_names = genob.iloc[:, 0].str.replace("\xa0", " ").str.replace("Â", "")
dup = sample_names.duplicated()
log(f"merged {genob.shape[0]} samples; dropped {int(dup.sum())} duplicate names")
genob = genob.loc[~dup, :]; sample_names = sample_names[~dup]
genob = genob.iloc[:, 1:]

genomat = np.multiply(np.subtract(genob.to_numpy(), 1), -2)     # {0,1,2}, NaN
maf = np.multiply(np.nanmean(genomat, axis=0), 0.5)
keepm = ~((maf == 0) | (maf > 0.05))                            # 0 < maf <= 0.05
mnf = marker_names[np.asarray(keepm)]; g = genomat[:, np.asarray(keepm)]
hetm = np.nanmean(g == 1, axis=0); g = g[:, hetm < 0.02]; mnf = mnf[hetm < 0.02]
hetl = np.nanmean(g == 1, axis=1); g = g[hetl < 0.02, :]; snf = sample_names[hetl < 0.02]
lines = snf.astype(str).tolist(); markers = mnf.tolist()
log(f"File_S10 filtered set: {len(lines)} lines x {len(markers)} markers")

# ---- assemble outputs ------------------------------------------------------
gr = pd.DataFrame(np.where(np.isnan(g), 3, g).astype(int), index=lines, columns=markers)
assert set(np.unique(gr.values)).issubset({0, 1, 2, 3})

md = pd.DataFrame({"marker": markers})
md["chrom"] = md["marker"].str.replace(r"_.+$", "", regex=True).str.replace("^S", "", regex=True).astype(int)
md["pos"] = md["marker"].str.replace(r"^S\d+_", "", regex=True).astype(int)

n_markers = len(markers)
r = 2 * 1500 / (100 * n_markers)                               # File S16 avg_r
mr = float((gr.values == 3).mean())
params = dict(PARAMS, r=r, mr=mr, n_markers=n_markers, n_lines=len(lines))

s18 = pd.read_csv(S18, index_col=0); s18.index = s18.index.astype(str)
s18a = s18.loc[lines, markers]                                 # exact same set (verified)

gr.to_csv(os.path.join(OUT, "geno_recoded.csv"))
md.to_csv(os.path.join(OUT, "markers.csv"), index=False)
pd.Series(lines).to_csv(os.path.join(OUT, "lines.csv"), index=False, header=False)
s18a.to_csv(os.path.join(OUT, "s18_aligned.csv"))
with open(os.path.join(OUT, "params.json"), "w") as fh:
    json.dump(params, fh, indent=2)
log(f"wrote geno_recoded {gr.shape}; avg_r={r:.3e}, mr={mr:.4f}")
log("done.")
