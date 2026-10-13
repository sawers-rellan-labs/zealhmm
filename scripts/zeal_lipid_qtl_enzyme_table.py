#!/usr/bin/env python3
# =============================================================================
# ZEAL lipid QTL enzyme-candidate table: one row per (significant QTL x
# glycerolipid candidate gene falling in its LOD-support CI), with the v5
# assembly locus id, the locus symbol, the enzyme/function, and the pathway(s).
#
# Symbol/name/function source: the MaizeGDB curated annotation
# data/Zm00001eb.1.fulldata.txt (authoritative; carries paper-assigned symbols
# like hpc1 that Gramene/MaizeMine lack). Columns (no header, tab): 2=gene_id,
# 11=symbol, 12=full_name, 13=description(function), 14=GO. Uncharacterized rows
# carry a LOC##### placeholder in col11 and "-" in col12 -> treated as no symbol.
# Pathways: from the v5 candidate build (data/zeal/lipid_candidates_v5.tsv).
#
# In : results/sim/zeal/rqtl/lipid/zeal_lipid_qtl_summary_annotated<sfx>.csv
#      data/zeal/lipid_candidates_v5.tsv, data/Zm00001eb.1.fulldata.txt
# Out: results/sim/zeal/rqtl/lipid/zeal_lipid_qtl_enzyme_candidates<sfx>.csv
# TAG env selects the scan: "" (default, log-scale main run) or e.g. "raw" -> "_raw" suffix.
# Run: python3 scripts/zeal_lipid_qtl_enzyme_table.py        # log-scale main run
#      TAG=raw python3 scripts/zeal_lipid_qtl_enzyme_table.py # raw-scale run
# =============================================================================
import os, re, csv, urllib.parse, urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAG = os.environ.get("TAG", "")
SFX = "" if TAG in ("", "taxon") else "_" + TAG
MM = "https://maizemine.rnet.missouri.edu/maizemine/service/query/results"
PLACEHOLDER = re.compile(r"^(LOC|Zm00001|cl\d)")


def is_symbol(s):
    return bool(s) and s != "-" and not PLACEHOLDER.match(s)


def mm_desc(genes):
    """MaizeMine v5 description fallback, only for genes MaizeGDB leaves blank."""
    if not genes:
        return {}
    xml = ('<query model="genomic" view="Gene.primaryIdentifier Gene.description">'
           '<constraint path="Gene.primaryIdentifier" op="ONE OF">'
           + "".join(f"<value>{g}</value>" for g in genes) + "</constraint></query>")
    url = MM + "?" + urllib.parse.urlencode({"query": xml, "format": "tab", "size": 5000})
    out = {}
    try:
        with urllib.request.urlopen(url, timeout=90) as r:
            for ln in r.readlines():
                p = ln.decode().rstrip("\n").split("\t")
                if len(p) >= 2:
                    out[p[0]] = p[1].strip('"').strip()
    except Exception:
        pass
    return out
ANN = os.path.join(ROOT, f"results/sim/zeal/rqtl/lipid/zeal_lipid_qtl_summary_annotated{SFX}.csv")
CAND = os.path.join(ROOT, "data/zeal/lipid_candidates_v5.tsv")
FULLDATA = os.path.join(ROOT, "data/Zm00001eb.1.fulldata.txt")
OUT = os.path.join(ROOT, f"results/sim/zeal/rqtl/lipid/zeal_lipid_qtl_enzyme_candidates{SFX}.csv")


def load_maizegdb():
    """gene_id -> (symbol, full_name, function). Placeholder LOC/'-' -> ''."""
    sym, name, func = {}, {}, {}
    with open(FULLDATA) as fh:
        for ln in fh:
            p = ln.rstrip("\n").split("\t")
            if len(p) < 13 or not p[1].startswith("Zm00001eb"):
                continue
            g = p[1]
            s = p[10].strip()
            sym[g] = s if is_symbol(s) else ""
            name[g] = "" if p[11].strip() in ("", "-") else p[11].strip()
            func[g] = p[12].strip()
    return sym, name, func


def main():
    pw = {}
    with open(CAND) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            pw[row["gene_v5"]] = row["pathway_names"]
    sym, name, func = load_maizegdb()

    qtls = list(csv.DictReader(open(ANN)))
    cand_genes = {g for q in qtls if q["candidate_genes"] for g in q["candidate_genes"].split(";")}
    # MaizeMine v5 description fallback only where MaizeGDB has no function AND no name
    need = [g for g in cand_genes if not func.get(g) and not name.get(g)]
    fallback = mm_desc(need)
    for g in need:
        if fallback.get(g):
            func[g] = fallback[g]
    rows = []
    for q in qtls:
        if not q["candidate_genes"]:
            continue
        for g in q["candidate_genes"].split(";"):
            rows.append({
                "trait": q["trait"], "class": q["class"], "chr": q["chr"],
                "peak_mb": round(float(q["peak_mb"]), 2), "lod": round(float(q["lod"]), 1),
                "ci_mb": f'{round(float(q["ci_left_mb"]),1)}-{round(float(q["ci_right_mb"]),1)}',
                "locus_id_v5": g, "symbol": sym.get(g, ""),
                "gene_label": sym.get(g) or g,  # symbol, else the assembly locus id
                "gene_name": name.get(g, ""),
                "enzyme_function": func.get(g, ""), "pathways": pw.get(g, ""),
            })
    rows.sort(key=lambda r: (-r["lod"], r["locus_id_v5"]))
    cols = ["trait", "class", "chr", "peak_mb", "lod", "ci_mb",
            "locus_id_v5", "symbol", "gene_label", "gene_name", "enzyme_function", "pathways"]
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    distinct = sorted({r["locus_id_v5"] for r in rows})
    nsym = len({r["locus_id_v5"] for r in rows if r["symbol"]})
    print(f"wrote {OUT}: {len(rows)} QTL x candidate rows, "
          f"{len(distinct)} distinct enzymes, {nsym} with a symbol")
    wsym = max((len(r["symbol"]) for r in rows), default=1) or 1
    for r in rows:
        print(f'{r["trait"]:<10} chr{r["chr"]:<2} {r["peak_mb"]:>7} Mb LOD{r["lod"]:>5}  '
              f'{r["locus_id_v5"]}  {r["symbol"]:<{wsym}}  {r["enzyme_function"] or r["gene_name"]}')


if __name__ == "__main__":
    main()
