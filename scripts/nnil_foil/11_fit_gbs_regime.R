# Step 2: fit the nNIL-GBS degradation regime from the REAL genotype calls.
#
# The released GBS carries no read depths, and Zhong et al. report no coverage
# number -- only "low coverage" (= TASSEL-GBS's 0.5-3x design band, Glaubitz 2014)
# and a missing-data rate. So coverage is a LATENT, pinned by the missing rate:
#
#     missing = pi_floor + (1 - pi_floor) * E[exp(-k * lambda)]
#
# From a single scalar mr only the PRODUCT k*lambda is identified; we take the
# naive Poisson split (k = 1, so lambda is the mean depth in reads) and fix the
# floor pi_floor = 0.01 (the best-covered line sits at ~0.9% missing, so the floor
# is provably <=0.01). The gamma spread of per-sample lambda is fit from the
# per-line missing distribution, and lambda_mean is solved via the gamma MGF so the
# gamma-AVERAGED missing equals the observed mr (correcting the Jensen gap that a
# plug-in lambda_mean would leave). het:hom is caller-dependent (TASSEL's GL-like
# call) so it is NOT used to fit -- it is reported only as a target to reproduce
# later through the caller.
#
#   Rscript scripts/nnil_foil/11_fit_gbs_regime.R
# Output: data/nnil_foil/nnil_gbs_regime.json  (lambda_mean, shape, k_decay,
#         pi_floor, error) + observed mr / state fractions for the record.

suppressMessages({
  library(data.table)
  library(here)
  library(jsonlite)
})

GENO <- here::here("data/nnil_equiv/geno_recoded.csv")
OUT <- here::here("data/nnil_foil/nnil_gbs_regime.json")
PI_FLOOR <- 0.01 # fixed floor (naive); capped by the min per-line missing
K <- 1.0 # naive Poisson split: P(>=1 read) = 1 - e^{-lambda}
ERROR <- 0.01 # Zhong's germ (homozygote genotyping error)

lg <- function(...) cat(sprintf(...), "\n")

# ---- observed missingness + state composition ------------------------------
x <- fread(GENO)
G <- as.matrix(x[, -1])
mr <- mean(G == 3)
per_line <- rowMeans(G == 3)
state <- c(REF = mean(G == 0), HET = mean(G == 1), DONORHOM = mean(G == 2), MISS = mr)
lg(
  "observed: mr = %.4f  | per-line min/mean/max = %.4f / %.4f / %.4f",
  mr, min(per_line), mean(per_line), max(per_line)
)
lg(
  "state fractions 0/1/2/3 = %.4f / %.4f / %.4f / %.4f  (het:hom = 1:%.1f, caller-dependent)",
  state["REF"], state["HET"], state["DONORHOM"], state["MISS"], state["DONORHOM"] / state["HET"]
)
if (PI_FLOOR > min(per_line)) {
  lg("NOTE: floor %.3f exceeds the min per-line missing %.4f (1 line); negligible.", PI_FLOOR, min(per_line))
}

# ---- per-line implied lambda -> gamma shape --------------------------------
# lambda_i solves m_i = pi + (1-pi) e^{-k lambda_i}; keep lines above the floor.
mi <- per_line[per_line > PI_FLOOR]
lam_i <- -log((mi - PI_FLOOR) / (1 - PI_FLOOR)) / K
shape <- mean(lam_i)^2 / var(lam_i) # method-of-moments gamma shape

# ---- solve lambda_mean so the gamma-averaged missing == mr -----------------
# E[e^{-k lambda}] = (1 + k*scale)^{-shape}, scale = lambda_mean/shape.
# target T = (mr - pi)/(1-pi):  lambda_mean = shape * (T^{-1/shape} - 1) / k
T <- (mr - PI_FLOOR) / (1 - PI_FLOOR)
lambda_mean <- shape * (T^(-1 / shape) - 1) / K
kl <- -log(T) # the identified product k*lambda (naive plug-in, for the record)
lg("fit: k = %.1f, pi_floor = %.2f  =>  k*lambda = %.3f", K, PI_FLOOR, kl)
lg("     gamma shape = %.3f, lambda_mean = %.3f reads/site (MGF-corrected)", shape, lambda_mean)

# ---- verify: Monte-Carlo the missing process, check it recovers mr ---------
set.seed(1)
N <- 200000L
lam <- pmax(0.01, rgamma(N, shape = shape, scale = lambda_mean / shape))
present <- (1 - PI_FLOOR) * (1 - exp(-K * lam))
sim_mr <- mean(1 - present)
lg("verify: simulated mean missing = %.4f  (target mr = %.4f; abs err %.4f)", sim_mr, mr, abs(sim_mr - mr))

regime <- list(
  source = "nnil_gbs", lambda_mean = lambda_mean, shape = shape,
  pi_floor = PI_FLOOR, k_decay = K, error = ERROR,
  observed_mr = mr, k_lambda = kl,
  target_state_fractions = as.list(round(state, 5)),
  note = "coverage inferred from mr (no depths released); het:hom is caller-dependent, not fit"
)
write_json(regime, OUT, auto_unbox = TRUE, pretty = TRUE, digits = 6)
lg("wrote %s", OUT)
