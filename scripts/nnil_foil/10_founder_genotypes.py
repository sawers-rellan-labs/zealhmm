#!/usr/bin/env python
# Step 1 (count-from-parents sim): extract the 24 NAM-founder chip genotypes on
# the SAME v5 GBS marker frame as 02_chip_truth.py's chip_truth_projected.csv, so
# the founders are co-registered with the chip truth and the real GBS data. This
# is the parent-genotype input to scripts/simulate_nested_nils.R, where donor
# non-informativeness (a founder carrying the B73/REF allele) emerges mechanically
# instead of via a scalar `nir` knob.
#
# Mapping mirrors 02_chip_truth.py VERBATIM (v3->v4 lift by name, chr 1-10, V4
# sort, recode 0/0,0/1,1/1,./. -> 0,1,2,3, same B73-REF + NIL-maf<0.2 QC, project
# to nearest GBS marker, relabel v4->v5). The ONLY difference: we keep the founder
# rows (not NIL) and do NOT run the HMM.
#
#   python scripts/nnil_foil/10_founder_genotypes.py
# Output (data/nnil_foil/):
#   founders_v5.csv   founder lines x v5-GBS-marker columns, genotype {0,1,2,3}
#                     (0=REF/B73, 1=het, 2=ALT-hom, 3=missing); REF=B73 by design

import os
import pandas as pd
import numpy as np

ROOT = "/Users/fvrodriguez/repos/zealhmm"
NNIL = os.path.join(ROOT, "agent/nNIL")
OUT = os.path.join(ROOT, "data/nnil_foil")


def log(msg):
    print(f"[10_founder_genotypes] {msg}", flush=True)


# ---- read chip genotypes (File_S02) ----------------------------------------
xlsx = os.path.join(NNIL, "File_S02.Chip data of NAM parents and nNILs v2.xlsx")
chip_header = pd.read_excel(xlsx, header=None, nrows=1)
chip_lines = chip_header.loc[0, 9:]
chip = pd.read_excel(xlsx, header=0, skiprows=1)
chip.rename(columns={"#CHROM": "chr", "POS (V3)": "pos_V3"}, inplace=True)
chip.rename(columns=dict(zip(chip.columns[9:], chip_lines)), inplace=True)
log(f"chip xlsx: {chip.shape[0]} SNPs x {len(chip_lines)} lines")

# ---- lift chip SNP positions V3 -> V4 (by shared marker name) --------------
v3 = pd.read_table(os.path.join(NNIL, "File_S12.nNIL_chip_SNP_positions_v3_6col.bed"), header=None)
v3.columns = ["chr", "startV3", "pos_V3", "name", "score", "strand"]
v3.drop(columns=["startV3", "score", "strand"], inplace=True)
v4 = pd.read_table(os.path.join(NNIL, "File_S13.nNIL_chip_SNP_positions_converted_to_V4.bed"), header=None)
v4.columns = ["chr_V4", "startV4", "pos_V4", "name", "score", "strand"]
v4.drop(columns=["startV4", "score", "strand"], inplace=True)
v3to4 = pd.merge(v3, v4, how="inner", on="name")
v3to4 = v3to4.loc[v3to4.chr_V4.isin([str(i) for i in range(1, 11)])]
v3to4["chr_V4"] = v3to4["chr_V4"].astype("int")

chip = chip.loc[chip.chr.isin([str(i) for i in range(1, 11)])]
chip = chip.astype({"chr": "int64"})
chip = chip.merge(v3to4, on=["chr", "pos_V3"])
chip.sort_values(by=["chr_V4", "pos_V4"], inplace=True)  # V4 order is critical
chip["name"] = "S" + chip["chr_V4"].astype("str") + "_" + chip["pos_V4"].astype("str")

# transpose: lines x markers
chipT = chip.iloc[:, 9 : (chip.shape[1] - 3)].transpose()
chipT.columns = chip["name"]
chipT.drop(
    index=["H100", "Ki3", "NC262", "NC304", "DRIL32.90 ", "DRIL32.095",
           "DRIL52.055", "DRIL62.078", "Mo17", "NIL-1030", "NIL-1290"],
    inplace=True,
)


def converter(x):
    if x == "0/0":
        return "0"
    if x == "0/1" or x == "1/0":
        return "1"
    if x == "1/1":
        return "2"
    if x == "./." or x == "nan" or pd.isna(x):
        return "3"


chipNp = chipT.applymap(converter).to_numpy(dtype="int64")

# ---- marker QC (File_S14 verbatim: same kept-marker set as chip_truth) ------
B73_nobs = np.sum(chipNp[0:2, ] != 3, axis=0)
B73_nobs[B73_nobs == 0] = 1
B73_afs = np.divide(
    np.sum(chipNp[0:2, ], axis=0) - 3 * np.sum(chipNp[0:2, ] == 3, axis=0), (2 * B73_nobs)
)
chipNIL = chipNp[chipT.index.str.contains(r"(NIL)"), ]
NIL_nobs = np.sum(chipNIL != 3, axis=0)
NIL_nobs[NIL_nobs == 0] = 1
NIL_afs = np.divide(
    np.sum(chipNIL, axis=0) - 3 * np.sum(chipNIL == 3, axis=0), (2 * NIL_nobs)
)
keep = np.logical_and(B73_afs == 0, NIL_afs < 0.20)
chipNp = chipNp[1::, keep]                       # drop one B73 (row 0), keep QC markers
chipSamples = chipT.index.to_series()[1:]
chipMarkers = chipT.columns.to_series()[keep]    # v4 names
log(f"after QC: {chipNp.shape[0]} lines x {chipNp.shape[1]} chip markers")

# ---- founders = NOT NIL and NOT B73 (the 24 NAM donors) ---------------------
is_founder = (~chipSamples.str.contains(r"(NIL)") & ~chipSamples.str.contains(r"B73")).to_numpy()
founders = chipNp[is_founder, :]
founder_names = chipSamples[is_founder].tolist()
log(f"genotyped NAM founders on chip: {founders.shape[0]} -> {founder_names}")

# ---- restrict to the ACTUAL pedigree donors of the 888-line GBS population ---
# The chip genotyped 24 NAM founders, but only 18 were used to create the nNILs
# (the extras are candidate donors for resolving pedigree mismatches). Keep only
# the pedigree donors, derived from the GBS NIL names ("<donor>/B73 NIL-xxxx").
lines = pd.read_csv(os.path.join(ROOT, "data/nnil_equiv/lines.csv"), header=None)
ped_donors = sorted(set(lines.iloc[:, 0].astype(str).str.replace(r"/B73.*", "", regex=True)))
keep_f = [i for i, n in enumerate(founder_names) if n in ped_donors]
missing = sorted(set(ped_donors) - set(founder_names))
founders = founders[keep_f, :]
founder_names = [founder_names[i] for i in keep_f]
log(f"pedigree donors in population: {len(ped_donors)} -> {ped_donors}")
log(f"kept {len(founder_names)} pedigree donors with chip genotypes; "
    f"MISSING (no chip genotype, get from HapMap): {missing}")

# ---- project onto nearest GBS marker + relabel v4->v5 (mirror 02) -----------
xw2g = pd.read_csv(os.path.join(OUT, "chip_markers_to_gbs.csv"))  # name(v4), closestGBS(v4)
chipMarkersDF = xw2g.loc[xw2g["name"].isin(chipMarkers)].copy()
uniq = chipMarkersDF.groupby(["chr", "closestGBS"], as_index=False).first()
uniq.sort_values(by=["chr", "pos"], inplace=True)

fdf = pd.DataFrame(founders, index=founder_names, columns=chipMarkers)
fdf = fdf.loc[:, fdf.columns.isin(uniq["name"])]
fdf.columns = uniq.set_index("name").loc[fdf.columns, "closestGBS"].values  # v4 GBS names

xwalk = pd.read_table(os.path.join(OUT, "markers_v5.tsv"))
v4_to_v5 = dict(zip(xwalk["marker_v4"], xwalk["marker"]))
mapped = fdf.columns.to_series().isin(v4_to_v5)
fdf = fdf.loc[:, mapped.values]
fdf.columns = [v4_to_v5[c] for c in fdf.columns]  # v5 ids

# align to the chip-truth frame (identical marker set as the evaluation target)
chip_truth = pd.read_csv(os.path.join(OUT, "chip_truth_projected.csv"), index_col=0)
shared = [c for c in chip_truth.columns if c in fdf.columns]
fdf = fdf.loc[:, shared]
fdf.index.name = "founder"
fdf.to_csv(os.path.join(OUT, "founders_v5.csv"))
log(f"wrote founders_v5.csv : {fdf.shape[0]} founders x {fdf.shape[1]} v5 markers "
    f"(aligned to chip_truth_projected frame)")

# ---- sanity: per-founder non-informative rate f0 (should be ~0.59) ----------
g = fdf.to_numpy()
f0s = [np.mean(row[row != 3] == 0) for row in g]
log(f"per-founder f0 (REF/non-informative): mean {np.mean(f0s):.3f}  "
    f"range [{np.min(f0s):.3f}, {np.max(f0s):.3f}]  (Holland grid nir=0.90)")