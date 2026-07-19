#!/usr/bin/env python
"""Timing worker: Holland File_S11 caller at one thinned marker size (odd-index
thinning, as in the RTIGER sweep), over the FULL 884-line population, from the
memory-mapped .bed. Peak RSS via parent /usr/bin/time -l.
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
val = open_bed(os.path.join(OUT, "geno.bed")).read()          # 884 x 64025, {0,1,2,NaN}

idx = np.arange(val.shape[1])                                  # odd-index thin
for _ in range(a.level):
    idx = idx[::2]
geno = np.where(np.isnan(val[:, idx]), 3, val[:, idx]).astype(float)
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
