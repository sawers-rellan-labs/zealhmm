#!/usr/bin/env python
"""
Run Holland's ORIGINAL nNIL caller (File_S11 call_intros, hmmlearn) on the shared
recoded genotype table from 01_build_input.py, with the published operating point.

  ~/anaconda3/envs/nilhmm/bin/python scripts/nnil_equiv/02_holland_calls.py

Output: data/nnil_equiv/holland_calls.csv  (lines x markers, states {0,1,2}).
"""
import importlib.util
import json
import os
import numpy as np
import pandas as pd
from bed_reader import open_bed

ROOT = os.path.expanduser("~/repos/zealhmm")
OUT = os.path.join(ROOT, "data/nnil_equiv")
S11 = os.path.join(ROOT, "agent/nNIL/File_S11_callIntrogressions.py")
BED = os.path.join(OUT, "geno.bed")


def log(m):
    print(f"[02_holland] {m}", flush=True)


# import Holland's caller module
spec = importlib.util.spec_from_file_location("cintro", S11)
ci = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ci)

# Stream genotypes from the compact .bed (same lean path as nilHMM's harness).
val = open_bed(BED).read()                               # lines x markers, {0,1,2,NaN}
geno = np.where(np.isnan(val), 3, val).astype(float)     # NaN -> 3
md = pd.read_csv(os.path.join(OUT, "markers.csv"))       # marker + chrom, .bed column order
params = json.load(open(os.path.join(OUT, "params.json")))
markers = md["marker"].tolist()
log(f"input(.bed): {geno.shape[0]} lines x {geno.shape[1]} markers; r={params['r']:.3e} mr={params['mr']:.4f}")

# marker_dict: chrom -> integer column indices (.bed column order == markers.csv)
marker_dict = {c: md.index[md["chrom"] == c].to_numpy() for c in range(1, 11)}

calls = ci.call_intros(
    geno=geno, marker_dict=marker_dict,
    nir=params["nir"], germ=params["germ"], gert=params["gert"], p=params["p"],
    mr=params["mr"], r=params["r"], f_1=params["f_1"], f_2=params["f_2"],
    return_calls=True,
)
lines = [l.strip() for l in open(os.path.join(OUT, "lines.csv"))]
pd.DataFrame(calls, index=lines, columns=markers).to_csv(os.path.join(OUT, "holland_calls.csv"))
vc = pd.Series(calls.ravel()).value_counts(normalize=True).round(4).to_dict()
log(f"wrote holland_calls.csv {calls.shape}; state fractions {vc}")
