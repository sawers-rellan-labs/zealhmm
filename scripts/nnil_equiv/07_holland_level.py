#!/usr/bin/env python
"""Equivalence worker: Holland's File_S11 caller at one thinned marker size
(odd-index thinning, as in the RTIGER sweep), over the FULL population, from the
memory-mapped .bed. Unlike 05_holland_worker.py (timing only) this one SAVES the
decoded calls so 07_equiv_sweep.R can compare them position-by-position against
nilHMM's nnil at the same size.

Writes the call matrix as a C-order int8 raw binary (lines x markers); R reads it
back with readBin. Both cores use the identical per-size r = 2*1500/(100*markers),
so any residual difference is a decoding difference, not a parameter difference.

Usage: python 07_holland_level.py --level L --out <path.bin>
"""
import argparse
import importlib.util
import json
import os
import numpy as np
import pandas as pd
from bed_reader import open_bed

ROOT = os.path.expanduser("~/repos/zealhmm")
OUT = os.path.join(ROOT, "data/nnil_equiv")
S11 = os.path.join(ROOT, "agent/nNIL/File_S11_callIntrogressions.py")

ap = argparse.ArgumentParser()
ap.add_argument("--level", type=int, default=0)
ap.add_argument("--out", required=True)
a = ap.parse_args()

spec = importlib.util.spec_from_file_location("cintro", S11)
ci = importlib.util.module_from_spec(spec); spec.loader.exec_module(ci)

DIR = os.path.join(OUT, f"thin_L{a.level}")                   # PRE-SPLIT: only this level
md = pd.read_csv(os.path.join(DIR, "markers.csv"))
params = json.load(open(os.path.join(OUT, "params.json")))
val = open_bed(os.path.join(DIR, "geno.bed")).read()          # lines x M (already thinned), {0,1,2,NaN}

geno = np.where(np.isnan(val), 3, val).astype(float)
chrom = md["chrom"].to_numpy()
marker_dict = {c: np.where(chrom == c)[0] for c in range(1, 11)}
r = 2 * 1500 / (100 * geno.shape[1])                          # per-size avg_r

calls = ci.call_intros(geno=geno, marker_dict=marker_dict, nir=params["nir"], germ=params["germ"],
                       gert=params["gert"], p=params["p"], mr=params["mr"], r=r,
                       f_1=params["f_1"], f_2=params["f_2"], return_calls=True)
np.ascontiguousarray(calls, dtype=np.int8).tofile(a.out)
print(f"RESULT caller=holland level={a.level} markers={geno.shape[1]} lines={geno.shape[0]} out={a.out}",
      flush=True)
