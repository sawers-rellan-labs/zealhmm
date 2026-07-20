#!/usr/bin/env python
"""Timing worker: Holland File_S11 caller at one thinned marker size (odd-index
thinning, as in the RTIGER sweep), over the full 888-line population, from the
memory-mapped .bed. Peak RSS via parent /usr/bin/time -l; reads ONLY the thinned
columns so peak RSS reflects the caller's working set at that density, not a
full-panel read.
Usage: python 05_holland_worker.py --level L"""
import argparse
import importlib.util
import json
import os
import time
import numpy as np
import pandas as pd
from bed_reader import open_bed

ROOT = os.path.expanduser("~/repos/zealhmm")
OUT = os.path.join(ROOT, "data/nnil_equiv")
S11 = os.path.join(ROOT, "agent/nNIL/File_S11_callIntrogressions.py")

ap = argparse.ArgumentParser()
ap.add_argument("--level", type=int, default=0)
a = ap.parse_args()

spec = importlib.util.spec_from_file_location("cintro", S11)
ci = importlib.util.module_from_spec(spec); spec.loader.exec_module(ci)

md = pd.read_csv(os.path.join(OUT, "markers.csv"))
params = json.load(open(os.path.join(OUT, "params.json")))

idx = np.arange(md.shape[0])                                   # odd-index thin (from marker count)
for _ in range(a.level):
    idx = idx[::2]
# read ONLY the thinned columns from disk -> peak RSS scales with the density
val = open_bed(os.path.join(OUT, "geno.bed")).read(index=np.s_[:, idx])  # N x len(idx), {0,1,2,NaN}
geno = np.where(np.isnan(val), 3, val).astype(float)
chrom = md["chrom"].to_numpy()[idx]
marker_dict = {c: np.where(chrom == c)[0] for c in range(1, 11)}
r = 2 * 1500 / (100 * len(idx))                               # per-size avg_r

t0 = time.perf_counter()
ci.call_intros(geno=geno, marker_dict=marker_dict, nir=params["nir"], germ=params["germ"],
               gert=params["gert"], p=params["p"], mr=params["mr"], r=r,
               f_1=params["f_1"], f_2=params["f_2"], return_calls=True)
dt = time.perf_counter() - t0
print(f"RESULT caller=holland level={a.level} markers={geno.shape[1]} lines={geno.shape[0]} "
      f"seconds={dt:.4f}", flush=True)
