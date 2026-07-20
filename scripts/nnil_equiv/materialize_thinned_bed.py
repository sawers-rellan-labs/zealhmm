#!/usr/bin/env python
"""SPLIT step for the nNIL sweeps: write an odd-index-thinned PLINK .bed (+ markers.csv)
per density level, so every sweep worker loads ONLY the dataset it benchmarks and does
no in-script thinning (mirrors scripts/materialize_thinned_panel.R for RTIGER).

  ~/anaconda3/envs/nilhmm/bin/python scripts/nnil_equiv/materialize_thinned_bed.py

Writes data/nnil_equiv/thin_L{0..6}/geno.{bed,bim,fam} + markers.csv (gitignored).
Level L keeps every marker after L odd-index halvings (L0 = full 64,025 -> L6 ~1,001).
"""
import os
import numpy as np
import pandas as pd
from bed_reader import open_bed, to_bed

ROOT = os.path.expanduser("~/repos/zealhmm")
OUT = os.path.join(ROOT, "data/nnil_equiv")
LEVELS = range(0, 7)

md = pd.read_csv(os.path.join(OUT, "markers.csv"))        # marker, chrom (.bed column order)
lines = [l.strip() for l in open(os.path.join(OUT, "lines.csv"))]
val = open_bed(os.path.join(OUT, "geno.bed")).read()      # N x 64025, {0,1,2,NaN}

for level in LEVELS:
    idx = np.arange(md.shape[0])
    for _ in range(level):
        idx = idx[::2]
    d = os.path.join(OUT, f"thin_L{level}")
    os.makedirs(d, exist_ok=True)
    props = {
        "iid": np.asarray(lines, dtype="U"),
        "sid": md["marker"].to_numpy()[idx].astype("U"),
        "chromosome": md["chrom"].to_numpy()[idx].astype("U"),
    }
    to_bed(os.path.join(d, "geno.bed"), val[:, idx], properties=props)
    md.iloc[idx].to_csv(os.path.join(d, "markers.csv"), index=False)
    print(f"[materialize_bed] L{level}: {len(idx)} markers -> {d}", flush=True)
