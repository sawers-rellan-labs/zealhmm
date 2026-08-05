#!/usr/bin/env bash
# Do the "non-panel-ALT" reads (the ones extract_gatk_counts_50k.sh zeroes) recur across samples?
#   - if each such site appears in ~1 sample, scattered across bases -> consistent with random error
#   - if sites appear in MANY samples with the SAME off-panel base   -> a real locus (a genuine third
#     allele the biallelic panel dropped, or a mapping artifact), NOT error
#
# One streaming pass over all 1,439 samples in the GATK table. For each 50K-set site we count:
#   cov_samples  = samples with any read (ref+alt>0)
#   np_samples   = samples where GATK's alt base is neither N nor the panel ALT, with ALT_COUNT>0
#   np_reads     = total such reads
#   modal_np_base + its fraction of np_samples = is the off-panel base consistent?
#
#   bash scripts/check_nonpanel_alt_recurrence.sh
#   env: BZEASEQ_DIR (default mount), PANEL (default data/zeal/snp50k_alleles.tsv)
#   out: agent/nonpanel_alt_recurrence.tsv (per-site) + a printed summary
set -uo pipefail
BZ=${BZEASEQ_DIR:-/Volumes/rsstu/users/r/rrellan/BZea/bzeaseq}
IN="$BZ/50K/allelic_counts50K.tsv"
PANEL=${PANEL:-data/zeal/snp50k_alleles.tsv}
OUTFILE=${OUTFILE:-agent/nonpanel_alt_recurrence.tsv}

[ -r "$PANEL" ] || { echo "ERROR: no $PANEL" >&2; exit 1; }
[ -r "$IN" ]    || { echo "ERROR: cannot read $IN (mount down?)" >&2; exit 1; }

awk -F'\t' -v OUTFILE="$OUTFILE" '
  NR==FNR { palt[$1 SUBSEP $2] = $4; next }           # pass 1: snp50k_alleles.tsv (chr pos ref alt)
  FNR==1  { next }                                     # pass 2: header of the big table
  {
    site = $2 SUBSEP $3
    if (!(site in palt)) next                          # off the 50K set
    if ($4 + $5 > 0) covsamp[site]++
    if ($7 != "N" && $7 != palt[site] && $5 + 0 > 0) { # a non-panel-ALT read event
      np[site]++; npreads[site] += $5
      bk = site SUBSEP $7; basecnt[bk]++
      if (basecnt[bk] > modcnt[site]) { modcnt[site] = basecnt[bk]; modbase[site] = $7 }
    }
  }
  END {
    print "chr\tpos\tpanel_alt\tcov_samples\tnp_samples\tnp_reads\tmodal_np_base\tmodal_frac" > OUTFILE
    for (s in np) {
      split(s, a, SUBSEP)
      frac = np[s] > 0 ? modcnt[s] / np[s] : 0
      printf "%s\t%s\t%s\t%d\t%d\t%d\t%s\t%.2f\n", a[1], a[2], palt[s], covsamp[s], np[s], npreads[s], modbase[s], frac >> OUTFILE
      n = np[s]
      if (n==1) h1++; else if (n==2) h2++; else if (n<=5) h5++; else if (n<=10) h10++; else if (n<=50) h50++; else if (n<=200) h200++; else hbig++
      tot_sites++; tot_ev += n; tot_rd += npreads[s]
    }
    printf "\n=== non-panel-ALT recurrence over all 1,439 samples (50K set) ===\n"
    printf "sites with >=1 non-panel-ALT read event: %d\n", tot_sites+0
    printf "total (sample,site) events: %d ; total non-panel reads: %d\n", tot_ev+0, tot_rd+0
    printf "\ndistribution of sites by # samples showing a non-panel-ALT read:\n"
    printf "  1 sample      : %d\n", h1+0
    printf "  2 samples     : %d\n", h2+0
    printf "  3-5 samples   : %d\n", h5+0
    printf "  6-10 samples  : %d\n", h10+0
    printf "  11-50 samples : %d\n", h50+0
    printf "  51-200 samples: %d\n", h200+0
    printf "  >200 samples  : %d\n", hbig+0
    printf "\n(error-like => mass at 1 sample, low modal_frac; real-locus => a tail in many samples with high modal_frac)\n"
  }
' "$PANEL" "$IN"

echo; echo "=== top 20 most-recurrent non-panel-ALT sites (np_samples desc) ==="
{ head -1 "$OUTFILE"; tail -n +2 "$OUTFILE" | sort -t$'\t' -k5,5nr | head -20; } | column -t -s$'\t'
