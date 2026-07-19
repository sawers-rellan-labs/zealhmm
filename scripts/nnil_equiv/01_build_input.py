#!/usr/bin/env python
"""
Build the faithful shared nNIL genotype input for the equivalence check.

The genotype table fed to BOTH callers is Holland's ACTUAL analysis input: File S1
raw GBS genotypes, subset to exactly the markers and lines that appear in the
PUBLISHED calls (File S18). File S18's 64,025-marker x 888-line set *is* the
definitive filtered set, so subsetting to it reproduces Holland's filtered input
without re-deriving his maf/het filters. Recode {0,0.5,1,NA} -> {0,1,2,3} exactly
as File_S10 / File_S16 (missing -> 3), and stamp the published operating point.

Run with the existing `nilhmm` conda env:
  ~/anaconda3/envs/nilhmm/bin/python scripts/nnil_equiv/01_build_input.py [--n-lines N]

Outputs (data/nnil_equiv/, gitignored):
  geno_recoded.csv  lines x markers, values {0,1,2,3}
  markers.csv       marker, chrom, pos  (in genomic order per chrom)
  s18_aligned.csv   File S18 published calls, same lines x markers order
  params.json       nir, germ, gert, p, mr, r, f_1, f_2
"""
import argparse
import json
import os
import sys
import numpy as np
import pandas as pd

ROOT = os.path.expanduser("~/repos/zealhmm")
ZH = os.path.join(ROOT, "data/zhong2025")
NNIL = os.path.join(ROOT, "agent/nNIL")
OUT = os.path.join(ROOT, "data/nnil_equiv")
os.makedirs(OUT, exist_ok=True)

S1 = os.path.join(ZH, "File_S01.nNIL_raw_SNPs_bgi_id_miss20.txt")
S18 = os.path.join(ZH, "File_S18.nNIL_gbs_HMM_introgressionCalls_full_set.csv")
S9 = os.path.join(NNIL, "File_S09.bgi_nil_id.txt")

# Published operating point (File_S16 final call); f_1/f_2 = expected het / homoz-intro
PARAMS = dict(nir=0.9, germ=0.01, gert=0.0001, p=0.9, f_1=0.007813, f_2=0.011179)


def log(msg):
    print(f"[01_build_input] {msg}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-lines", type=int, default=0,
                    help="subset to first N common lines (0 = all)")
    args = ap.parse_args()

    # --- raw genotypes: skip <Numeric>, header row is <Marker> + marker names ---
    log("reading File S1 (raw genotypes)...")
    geno = pd.read_table(S1, sep="\t", skiprows=1, header="infer")
    geno = geno.rename(columns={geno.columns[0]: "bgi_id"}).set_index("bgi_id")
    log(f"  raw: {geno.shape[0]} samples x {geno.shape[1]} markers")

    # --- bgi_id -> nNIL line-name translation (File S9) ---
    trans = pd.read_table(S9, header=0, encoding="windows-1252")
    log(f"  File S9 columns: {list(trans.columns)}")
    bcol = [c for c in trans.columns if "bgi" in c.lower()][0]
    ncol = [c for c in trans.columns if c != bcol][0]
    bgi2name = dict(zip(trans[bcol].astype(str), trans[ncol].astype(str)))

    # --- published calls (File S18): the definitive filtered marker/line set ---
    log("reading File S18 (published calls)...")
    s18 = pd.read_csv(S18, index_col=0)
    s18.index = s18.index.astype(str)
    log(f"  S18: {s18.shape[0]} lines x {s18.shape[1]} markers")

    # translate geno index to line names; keep those present in S18
    geno.index = [bgi2name.get(str(b), str(b)) for b in geno.index]
    geno = geno[~geno.index.duplicated(keep="first")]      # drop dup translations

    common_lines = [ln for ln in s18.index if ln in set(geno.index)]
    common_markers = [m for m in s18.columns if m in set(geno.columns)]
    log(f"  common: {len(common_lines)} lines, {len(common_markers)} markers")
    if args.n_lines:
        common_lines = common_lines[: args.n_lines]
        log(f"  subset to first {len(common_lines)} lines")

    # marker order = genomic (chrom, pos) parsed from name S{chr}_{pos}
    md = pd.DataFrame({"marker": common_markers})
    md["chrom"] = md["marker"].str.replace(r"_.+$", "", regex=True).str.replace("^S", "", regex=True).astype(int)
    md["pos"] = md["marker"].str.replace(r"^S\d+_", "", regex=True).astype(int)
    md = md.sort_values(["chrom", "pos"]).reset_index(drop=True)
    markers = md["marker"].tolist()

    g = geno.loc[common_lines, markers].copy()
    s18a = s18.loc[common_lines, markers].copy()

    # recode {0,0.5,1,NA} -> {2,1,0,3}: (1-raw)*2, NA->3   (== File_S10/File_S16)
    gr = np.where(g.isna().values, 3, ((1.0 - g.values) * 2)).astype(int)
    gr = pd.DataFrame(gr, index=common_lines, columns=markers)
    assert set(np.unique(gr.values)).issubset({0, 1, 2, 3}), "recode produced unexpected values"

    # params derived exactly as File_S16
    n_markers = len(markers)
    r = 2 * 1500 / (100 * n_markers)                 # avg_r
    mr = float((gr.values == 3).mean())              # missing rate
    params = dict(PARAMS, r=r, mr=mr, n_markers=n_markers, n_lines=len(common_lines))

    gr.to_csv(os.path.join(OUT, "geno_recoded.csv"))
    md.to_csv(os.path.join(OUT, "markers.csv"), index=False)
    s18a.to_csv(os.path.join(OUT, "s18_aligned.csv"))
    with open(os.path.join(OUT, "params.json"), "w") as fh:
        json.dump(params, fh, indent=2)
    log(f"wrote geno_recoded {gr.shape}, markers {len(markers)}, avg_r={r:.3e}, mr={mr:.4f}")
    log(f"params: {params}")
    log("done.")


if __name__ == "__main__":
    main()