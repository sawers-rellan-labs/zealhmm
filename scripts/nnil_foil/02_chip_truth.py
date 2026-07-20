#!/usr/bin/env python
# Calibration foil, step 3: reproduce Holland's SNP-chip introgression calls (the
# gold standard for the foil) by porting File_S14 (Step 2) VERBATIM in logic and
# calling his own HMM (File_S11 ci.call_intros). The ONLY changes vs File_S14 are
# I/O paths and the GBS marker list source (our data/nnil_equiv/markers.csv, the
# exact File_S10-filtered 64,025-marker set). Nothing about the model is altered,
# so this is a reproduction, not an approximation.
#
# Pipeline (Holland File_S14): read chip xlsx -> lift chip SNP positions V3->V4
# (File_S12/S13, by shared marker name) -> chr 1-10, sort by V4 -> transpose to
# lines x markers -> drop the 11 non-analysis lines -> recode 0/0,0/1,1/1,./. to
# 0,1,2,3 -> QC (drop markers non-REF on B73 or NIL maf>0.2; drop one B73) -> chip
# HMM with Holland's chip-optimal params (nir=0.9 germ=0.01 gert=0.001 p=0.9
# r=avg_r/2) -> project calls onto nearest GBS marker names.
#
#   python scripts/nnil_foil/02_chip_truth.py
# Output (data/nnil_foil/):
#   chip_truth_projected.csv  NIL lines x GBS-marker columns, states {0,1,2}
#   chip_markers_to_gbs.csv   chip marker -> nearest GBS marker + distance

import os
import sys
import importlib.util
import pandas as pd
import numpy as np

ROOT = "/Users/fvrodriguez/repos/zealhmm"
NNIL = os.path.join(ROOT, "agent/nNIL")
OUT = os.path.join(ROOT, "data/nnil_foil")
os.makedirs(OUT, exist_ok=True)


def log(msg):
    print(f"[02_chip_truth] {msg}", flush=True)


# import Holland's HMM caller (File_S11) as the module `ci`
spec = importlib.util.spec_from_file_location(
    "ci", os.path.join(NNIL, "File_S11_callIntrogressions.py")
)
ci = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ci)

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

# drop the 11 lines Holland excluded from analysis (File_S14 verbatim)
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


chipT2 = chipT.applymap(converter)
chipNp = chipT2.to_numpy(dtype="int64")

# ---- marker QC (File_S14 verbatim) -----------------------------------------
# drop markers non-REF on the two B73 samples, and NIL maf > 0.2; drop one B73
B73_nobs = np.sum(chipNp[0:2,] != 3, axis=0)
B73_nobs[B73_nobs == 0] = 1
B73_afs = np.divide(
    np.sum(chipNp[0:2,], axis=0) - 3 * np.sum(chipNp[0:2,] == 3, axis=0), (2 * B73_nobs)
)
chipNIL = chipNp[chipT.index.str.contains(r"(NIL)"),]
NIL_nobs = np.sum(chipNIL != 3, axis=0)
NIL_nobs[NIL_nobs == 0] = 1
NIL_afs = np.divide(
    np.sum(chipNIL, axis=0) - 3 * np.sum(chipNIL == 3, axis=0), (2 * NIL_nobs)
)
keep = np.logical_and(B73_afs == 0, NIL_afs < 0.20)
chipNp = chipNp[1::, keep]
chipSamples = chipT.index.to_series()[1:]
chipMarkers = chipT.columns.to_series()[keep]
log(f"after QC: {chipNp.shape[0]} lines x {chipNp.shape[1]} chip markers")

# ---- nearest GBS marker to each chip marker (our marker list) --------------
gbsMarkers = pd.read_csv(os.path.join(ROOT, "data/nnil_equiv/markers.csv"))["marker"]
chrnames = gbsMarkers.str.extract(r"(S\d*)").iloc[:, 0].str.replace("S", "").astype("int64")
posnames = gbsMarkers.str.extract(r"(\d*$)").iloc[:, 0].astype("int64")
gbsMarkersDF = pd.DataFrame({"name": gbsMarkers, "chr": chrnames, "pos": posnames})

chipMarkersDF = chip[["name", "chr_V4", "pos_V4"]].copy()
chipMarkersDF = chipMarkersDF.loc[chipMarkersDF["name"].isin(chipMarkers)]
chipMarkersDF.rename(columns={"chr_V4": "chr", "pos_V4": "pos"}, inplace=True)
chipMarkersDF.chr = chipMarkersDF.chr.astype("int64")


def closestMatch(row):
    g = gbsMarkersDF.loc[gbsMarkersDF.chr == row["chr"], :]
    dist = abs(g.pos - row["pos"]).reset_index(drop=True)
    return g.iloc[dist.idxmin(), 0]


chipMarkersDF["closestGBS"] = chipMarkersDF.apply(closestMatch, axis=1)
chipMarkersDF["closest_gbs_pos"] = chipMarkersDF.closestGBS.str.extract(r"(\d*$)").astype("int64")
chipMarkersDF["distance"] = abs(chipMarkersDF["closest_gbs_pos"] - chipMarkersDF["pos"])
chipMarkersDF.to_csv(os.path.join(OUT, "chip_markers_to_gbs.csv"), index=False)
log(f"nearest-GBS distance median {chipMarkersDF.distance.median():.0f} bp, "
    f"max {chipMarkersDF.distance.max():.0f} bp")

# ---- chip HMM (Holland's chip-optimal parameters) --------------------------
donorCalls = chipNp[~chipSamples.str.contains(r"(NIL)") & ~chipSamples.str.contains(r"B73"), :].astype("float")
donorCalls[donorCalls == 3] = np.nan
maf = np.nanmean(np.multiply(np.nanmean(donorCalls, axis=1), 0.5))
nir = 1 - maf

chipNp = np.nan_to_num(chipNp, copy=True, nan=3)
missing_rate = (chipNp == 3).sum() / (chipNp.shape[0] * chipNp.shape[1])
avg_r = 2 * 1500 / (100 * len(chipMarkers))

chroms = chipMarkers.replace("_.+$", "", regex=True).replace("^S", "", regex=True).reset_index(drop=True)
markers_by_chrom = {i: chroms[chroms == str(i)].index for i in range(1, 11)}

f_2 = 0.011179
f_1 = 0.007813
# File_S14 chip-optimal: nir=0.9 germ=0.01 gert=0.001 p=0.9 r=avg_r/2
finalModel = ci.call_intros(
    geno=chipNp, marker_dict=markers_by_chrom,
    nir=0.9, germ=0.01, gert=0.001, p=0.9, mr=missing_rate, r=avg_r / 2,
    f_1=f_1, f_2=f_2, return_calls=True,
)
finalModel = pd.DataFrame(finalModel, index=chipSamples, columns=chipMarkers)
finalModel.index.name = "Line"

# ---- project NIL calls onto nearest GBS markers ----------------------------
NILIndices = chipSamples[chipSamples.str.contains(r"(NIL)") | (chipSamples == "B73 ")]
chipMarkersUnique = chipMarkersDF.groupby(["chr", "closestGBS"], as_index=False).first()
chipMarkersUnique.sort_values(by=["chr", "pos"], inplace=True)

proj = finalModel.loc[chipSamples.isin(NILIndices), chipMarkers.isin(chipMarkersUnique["name"])]
proj.columns = chipMarkersUnique["closestGBS"].values  # v4 GBS marker names
proj.index.name = "Line"

# Drop chip calls landing on GBS markers that did NOT lift to v5 (unmapped are
# absent from the crosswalk), and relabel the survivors to their v5 id.
xwalk = pd.read_table(os.path.join(OUT, "markers_v5.tsv"))  # marker(v5), marker_v4, ...
v4_to_v5 = dict(zip(xwalk["marker_v4"], xwalk["marker"]))
n_before = proj.shape[1]
mapped = proj.columns.to_series().isin(v4_to_v5)
proj = proj.loc[:, mapped.values]
proj.columns = [v4_to_v5[c] for c in proj.columns]  # v4 GBS name -> v5 id
log(f"projection: dropped {n_before - proj.shape[1]} unmapped GBS markers "
    f"({proj.shape[1]} kept, v5 ids)")

proj.to_csv(os.path.join(OUT, "chip_truth_projected.csv"))
log(f"wrote chip_truth_projected.csv : {proj.shape[0]} NIL lines x {proj.shape[1]} GBS-marker cols")
fr = pd.Series(proj.to_numpy().ravel()).value_counts(normalize=True).sort_index()
log(f"chip-truth state fractions: {dict(fr.round(4))}")
