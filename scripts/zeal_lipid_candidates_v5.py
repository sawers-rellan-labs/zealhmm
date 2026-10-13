#!/usr/bin/env python3
# =============================================================================
# ZEAL lipid candidate genes — B73 v5 rebuild of the glycerolipid list.
#
# Faithful reproduction of the pglipid method (sawers-rellan-labs/pglipid
# scripts/1_get_gplipid_genes.R) on the NEWEST versions, NOT a lift-over of the
# AGPv4 list. Controlled-vocabulary seed = KEGG maize pathways zma00561
# (Glycerolipid metabolism) + zma00564 (Glycerophospholipid metabolism); pathway
# expansion/curation over CornCyc, which on v5 is the MaizeMine "MaizeGDB E2P2
# pathways data set" (MetaCyc PWY- ids). All sources are Cloudflare-free:
#   - KEGG REST (live)                    : pathway -> Entrez gene ids
#   - MaizeMine InterMine REST (rel 1.6)  : Entrez -> v5 Zm00001eb (crossReferences)
#                                           v5 gene -> E2P2 pathways (Gene.pathways)
# (MaizeGDB www+ftp are Cloudflare-blocked from CLI; MaizeMine is the worked-out
# substitute — see memory maize-gene-id-mapping-maizemine.)
#
# pglipid logic reproduced:
#   seeds       = KEGG pathway genes mapped to v5
#   expand_cyc  = per KEGG pathway, union(seeds, genes of E2P2 pathways hit by >1 seed)
#   candidates  = union of expand over both KEGG pathways
#   curate      = E2P2 pathways fully covered by candidates (cover==1),
#                 minus REMOVE list, plus ADD (PWY-6825 = phosphatidylcholine biosyn V)
#   final       = union(all KEGG seed v5 genes, genes of curated pathways)
# (orphan-enzyme + KEGG-only-unannotated add-ons from the paper are NOT in the
#  script's written output there either; noted, not included here.)
#
# Out: data/zeal/lipid_candidates_v5.tsv  (gene_v5, chr, start, end, in_kegg_seed,
#        in_curated_pathway, n_glyco_pathways, pathway_ids, pathway_names)
#      data/zeal/lipid_candidates_v5_provenance.txt
# Run: python3 scripts/zeal_lipid_candidates_v5.py
# =============================================================================
import os, time, urllib.parse, urllib.request, collections

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
COORDS = os.path.join(ROOT, "data/ref/b73_gene_coords.tsv")
OUT = os.path.join(ROOT, "data/zeal/lipid_candidates_v5.tsv")
PROV = os.path.join(ROOT, "data/zeal/lipid_candidates_v5_provenance.txt")

KEGG_PATHWAYS = ["zma00561", "zma00564"]  # Glycerolipid, Glycerophospholipid
E2P2 = "MaizeGDB E2P2 pathways data set"
# manual curation, verbatim from pglipid 1_get_gplipid_genes.R
REMOVE = {"PWY-6446", "PWY-5143", "PWY66-21", "PWY-6457", "PWY-2501",
          "PWY-6333", "PWY1F-467", "PWY0-1313", "PWY-3941"}
ADD = {"PWY-6825"}  # phosphatidylcholine biosynthesis V

MM = "https://maizemine.rnet.missouri.edu/maizemine/service/query/results"
KEGG = "https://rest.kegg.jp"
log_lines = []


def log(msg):
    print(msg)
    log_lines.append(msg)


def http_lines(url, timeout=120):
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as r:
                return [ln.decode("utf-8", "replace").rstrip("\n") for ln in r.readlines()]
        except Exception:
            if attempt == 3:
                raise
            time.sleep(3 * (attempt + 1))


def kegg_entrez(pathway):
    # rest.kegg.jp/link/zma/<pathway> -> "path:zmaXXXXX\tzma:<entrez>"
    out = set()
    for ln in http_lines(f"{KEGG}/link/zma/{pathway}"):
        if "\t" in ln:
            g = ln.split("\t")[1]
            out.add(g.split(":", 1)[1] if ":" in g else g)
    return out


def mm_query(view, constraints, page=10000):
    """Run a paged MaizeMine query; return list of row (list-of-cols)."""
    cx = "".join(constraints)
    rows, start = [], 0
    while True:
        xml = f'<query model="genomic" view="{view}">{cx}</query>'
        url = MM + "?" + urllib.parse.urlencode(
            {"query": xml, "format": "tab", "start": start, "size": page})
        chunk = [ln.split("\t") for ln in http_lines(url) if ln]
        rows.extend(chunk)
        if len(chunk) < page:
            break
        start += page
    return rows


def strip_tx(g):
    # gene-level id; drop any transcript/protein suffix just in case
    for sep in ("_T", "_P"):
        if sep in g:
            g = g.split(sep)[0]
    return g


def main():
    # ---- 1. KEGG controlled-vocabulary seeds (live) --------------------------
    kegg_seed_entrez = {}
    for pw in KEGG_PATHWAYS:
        e = kegg_entrez(pw)
        kegg_seed_entrez[pw] = e
        log(f"KEGG {pw}: {len(e)} Entrez genes")
    all_entrez = set().union(*kegg_seed_entrez.values())
    log(f"KEGG union Entrez: {len(all_entrez)}")

    # ---- 2. Entrez -> v5 (MaizeMine crossReferences) -------------------------
    rows = mm_query(
        "Gene.primaryIdentifier Gene.crossReferences.identifier",
        ['<constraint path="Gene.organism.name" op="=" value="Zea mays"/>',
         '<constraint path="Gene.source" op="=" value="RefSeq"/>',
         '<constraint path="Gene.crossReferences.source.name" op="=" value="MaizeMine"/>'])
    entrez2v5 = collections.defaultdict(set)
    for r in rows:
        if len(r) >= 2 and r[1].startswith("Zm00001eb"):
            entrez2v5[r[0]].add(strip_tx(r[1]))
    log(f"MaizeMine RefSeq Entrez->v5 rows: {len(rows)}; mapped Entrez: {len(entrez2v5)}")
    mapped = {e for e in all_entrez if e in entrez2v5}
    log(f"KEGG Entrez with a v5 model: {len(mapped)}/{len(all_entrez)}")

    def to_v5(entrez_set):
        out = set()
        for e in entrez_set:
            out |= entrez2v5.get(e, set())
        return out

    kegg_seed_v5 = {pw: to_v5(kegg_seed_entrez[pw]) for pw in KEGG_PATHWAYS}
    all_kegg_seed_v5 = set().union(*kegg_seed_v5.values())
    log(f"KEGG seed v5 genes: {' '.join(f'{pw}={len(kegg_seed_v5[pw])}' for pw in KEGG_PATHWAYS)} union={len(all_kegg_seed_v5)}")

    # ---- 3. Full v5 gene <-> E2P2 (CornCyc) pathway table --------------------
    # E2P2 genes are keyed by primaryIdentifier = v5 Zm00001eb (source "B73 Zm00001eb.1"),
    # the SAME key the KEGG seeds were mapped into via crossReferences.
    rows = mm_query(
        "Gene.primaryIdentifier Gene.pathways.identifier Gene.pathways.name",
        ['<constraint path="Gene.organism.name" op="=" value="Zea mays"/>',
         f'<constraint path="Gene.pathways.dataSets.name" op="=" value="{E2P2}"/>'])
    pwy_genes = collections.defaultdict(set)
    gene_pwys = collections.defaultdict(set)
    pwy_name = {}
    for r in rows:
        if len(r) >= 2 and r[0].startswith("Zm00001eb") and r[1]:
            g, pid = strip_tx(r[0]), r[1]
            pwy_genes[pid].add(g)
            gene_pwys[g].add(pid)
            if len(r) >= 3:
                pwy_name[pid] = r[2]
    log(f"E2P2: {len(pwy_genes)} pathways, {len(gene_pwys)} v5 genes, {len(rows)} rows")

    # ---- 4. expand_cyc per KEGG pathway (pathways hit by >1 seed) ------------
    candidates = set()
    for pw in KEGG_PATHWAYS:
        seeds = kegg_seed_v5[pw]
        ntest = collections.Counter()
        for g in seeds:
            for pid in gene_pwys.get(g, ()):
                ntest[pid] += 1
        expand = set(seeds)
        expanded_pwys = [pid for pid, n in ntest.items() if n > 1]
        for pid in expanded_pwys:
            expand |= pwy_genes[pid]
        candidates |= expand
        log(f"expand {pw}: seeds={len(seeds)} pathways(n_test>1)={len(expanded_pwys)} -> expand={len(expand)}")
    log(f"candidates (union expand): {len(candidates)}")

    # ---- 5. curate: E2P2 pathways fully covered by candidates ----------------
    fully = [pid for pid, gs in pwy_genes.items() if gs and gs <= candidates]  # cover==1
    kept = (set(fully) - REMOVE) | (ADD & set(pwy_genes))
    missing_add = ADD - set(pwy_genes)
    log(f"pathways cover==1: {len(fully)}; after remove({len(REMOVE)})+add({len(ADD)}): {len(kept)}")
    if missing_add:
        log(f"  NOTE: ADD pathway(s) absent from v5 E2P2, cannot add: {sorted(missing_add)}")
    absent_remove = REMOVE - set(pwy_genes)
    if absent_remove:
        log(f"  NOTE: REMOVE ids absent from v5 E2P2 (no-op): {sorted(absent_remove)}")
    curated_genes = set()
    for pid in kept:
        curated_genes |= pwy_genes.get(pid, set())
    log(f"curated pathway genes: {len(curated_genes)} from {len(kept)} pathways")

    # ---- 6. final = union(KEGG seed v5, curated) ----------------------------
    final = all_kegg_seed_v5 | curated_genes
    log(f"FINAL v5 candidates: {len(final)}  (KEGG-seed {len(all_kegg_seed_v5)}, curated {len(curated_genes)}, overlap {len(all_kegg_seed_v5 & curated_genes)})")

    # ---- 7. annotate (coords + glycerolipid pathway membership) --------------
    coords = {}
    with open(COORDS) as fh:
        for ln in fh:
            p = ln.rstrip("\n").split("\t")
            if len(p) >= 4:
                coords[p[0]] = (p[1], p[2], p[3])
    kept_names = {pid: pwy_name.get(pid, "") for pid in kept}
    with open(OUT, "w") as fh:
        fh.write("gene_v5\tchr\tstart\tend\tin_kegg_seed\tin_curated_pathway\t"
                 "n_glyco_pathways\tpathway_ids\tpathway_names\n")
        for g in sorted(final):
            chrm, s, e = coords.get(g, ("", "", ""))
            gp = sorted(gene_pwys.get(g, set()) & kept)  # its curated glycerolipid pathways
            fh.write("\t".join([
                g, chrm, s, e,
                str(g in all_kegg_seed_v5), str(g in curated_genes),
                str(len(gp)), ";".join(gp),
                ";".join(kept_names.get(p, "") for p in gp)]) + "\n")
    n_nocoord = sum(1 for g in final if g not in coords)
    log(f"wrote {OUT} ({len(final)} genes; {n_nocoord} without local coords)")
    with open(PROV, "w") as fh:
        fh.write("\n".join(log_lines) + "\n")
    log(f"wrote {PROV}")


if __name__ == "__main__":
    main()
