#!/usr/bin/env python
"""Timing worker: Holland File_S11 caller over the full 888-line population, from the
PRE-SPLIT memory-mapped .bed for one density level (thin_L<level>/, written by
materialize_thinned_bed.py). Loads ONLY the level it benchmarks; no in-script
thinning, so peak RSS (parent /usr/bin/time -l) reflects the caller's working set at
that density.
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

DIR = os.path.join(OUT, f"thin_L{a.level}")                   # PRE-SPLIT: only this level
md = pd.read_csv(os.path.join(DIR, "markers.csv"))
params = json.load(open(os.path.join(OUT, "params.json")))

val = open_bed(os.path.join(DIR, "geno.bed")).read()          # N x M (already thinned), {0,1,2,NaN}
geno = np.where(np.isnan(val), 3, val).astype(float)
chrom = md["chrom"].to_numpy()
marker_dict = {c: np.where(chrom == c)[0] for c in range(1, 11)}
r = 2 * 1500 / (100 * geno.shape[1])                          # per-size avg_r

t0 = time.perf_counter()
ci.call_intros(geno=geno, marker_dict=marker_dict, nir=params["nir"], germ=params["germ"],
               gert=params["gert"], p=params["p"], mr=params["mr"], r=r,
               f_1=params["f_1"], f_2=params["f_2"], return_calls=True)
dt = time.perf_counter() - t0
print(f"RESULT caller=holland level={a.level} markers={geno.shape[1]} lines={geno.shape[0]} "
      f"seconds={dt:.4f}", flush=True)
