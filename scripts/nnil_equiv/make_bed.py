#!/usr/bin/env python
"""
Convert the recoded genotype table to PLINK .bed (compact, memory-mappable) so the
callers can stream rows instead of slurping a 114 MB wide CSV. Our encoding is
stored directly as the .bed value: 0 = REF-homoz (B73), 1 = het, 2 = donor-homoz,
missing = NaN. Both bed_reader (Python) and BEDMatrix (R) return these verbatim
(verified: no allele-orientation flip), so `0 = REF` is enforced here by a
ROUND-TRIP ASSERT -- write, read back, require identity -- not by trusting a
convention. Any orientation/missing-code slip fails loudly at staging.

  ~/anaconda3/envs/nilhmm/bin/python scripts/nnil_equiv/make_bed.py

Output: data/nnil_equiv/geno.{bed,bim,fam}
"""
import os
import numpy as np
import pandas as pd
from bed_reader import to_bed, open_bed

OUT = os.path.expanduser("~/repos/zealhmm/data/nnil_equiv")
BED = os.path.join(OUT, "geno.bed")


def log(m):
    print(f"[make_bed] {m}", flush=True)


gr = pd.read_csv(os.path.join(OUT, "geno_recoded.csv"), index_col=0)
md = pd.read_csv(os.path.join(OUT, "markers.csv")).set_index("marker").loc[gr.columns].reset_index()
lines, markers = gr.index.astype(str).tolist(), gr.columns.tolist()
log(f"recoded: {len(lines)} lines x {len(markers)} markers")

val = gr.to_numpy().astype(float)          # lines x markers, {0,1,2,3}
orig = val.copy()
val[val == 3] = np.nan                     # 3 -> missing for .bed

# Line names contain spaces (".../B73 NIL-1001"), and PLINK .fam is whitespace-
# delimited -> the IID would be mangled. Keep the REAL names in a sidecar (row
# order) and give the .fam a safe positional IID; callers align by row order.
pd.Series(lines).to_csv(os.path.join(OUT, "lines.csv"), index=False, header=False)
safe_iid = [f"L{i:04d}" for i in range(len(lines))]
to_bed(
    BED, val,
    properties={
        "iid": safe_iid, "sid": markers,
        "chromosome": md["chrom"].astype(str).tolist(),
        "bp_position": md["pos"].astype(int).tolist(),
        "allele_1": ["donor"] * len(markers),  # count-of-A1 == donor dose -> 0 = REF
        "allele_2": ["REF"] * len(markers),
    },
)
sz = os.path.getsize(BED) / 1024**2
log(f"wrote geno.bed ({sz:.1f} MiB) + .bim/.fam")

# ---- ROUND-TRIP ASSERT: read back, require identity to the original {0,1,2,3} ----
rb = open_bed(BED).read()                  # {0,1,2,NaN}
rb_recoded = np.where(np.isnan(rb), 3, rb).astype(int)
if not np.array_equal(rb_recoded, orig.astype(int)):
    n = int((rb_recoded != orig.astype(int)).sum())
    raise SystemExit(f"[make_bed] ROUND-TRIP FAILED: {n} cells differ after .bed write/read "
                     "-- encoding or allele orientation is wrong; refusing to proceed.")
log("round-trip OK: .bed reads back identical to the recoded table (0=REF preserved).")
