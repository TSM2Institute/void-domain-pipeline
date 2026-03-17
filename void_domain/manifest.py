# void_domain/manifest.py
# VOID-DOMAIN v1 — Sealed Manifest
# Locked: 15 March 2026
# DO NOT MODIFY after first commit without creating manifest_v2.py

# ── Void Catalogue ──────────────────────────────────────────────────────────
VOID_CATALOGUE_REF     = "Mao et al. (2017) ApJ 835 161 — SDSS DR12 BOSS"
VOID_N_EXPECTED        = 1228          # after quality cuts
VOID_REFF_MIN_HMPC     = 20.0          # h^-1 Mpc
VOID_REFF_MAX_HMPC     = 100.0         # h^-1 Mpc
VOID_Z_MIN             = 0.2
VOID_Z_MAX             = 0.7

# ── Cluster Catalogue ────────────────────────────────────────────────────────
CLUSTER_CATALOGUE_REF  = "redMaPPer DR8 — Rykoff et al. (2014) ApJ 785 104"
CLUSTER_MASS_PROXY     = "M200c from richness-lambda conversion"

# ── M_surrounding Definition (Geoffrey Thwaites, updated 17 March 2026) ──────
# PRIMARY: annulus method
ANNULUS_INNER_FACTOR   = 2.0           # inner edge = 2 × R_void
ANNULUS_OUTER_FACTOR   = 4.0           # outer edge = 4 × R_void
MASS_UNIT              = 1e14          # solar masses

M_SURROUNDING_DEFINITION = (
    "Density contrast delta = (rho_annulus - rho_mean) / rho_mean. "
    "rho_annulus = sum(M_200c in annulus) / V_annulus. "
    "Annulus: 2 x R_void to 4 x R_void from void centre. "
    "rho_mean = total group mass / box volume. "
    "Locked per Geoffrey Thwaites directive 17 March 2026."
)

# SECONDARY (parallel, robustness only)
SPHERE_FACTOR          = 3.0           # d < 3 × R_void sphere
# TERTIARY: density-shell integrated mass (implemented separately)

# ── R_void Definition ────────────────────────────────────────────────────────
# Use catalogue-published R_eff = (3V/4π)^(1/3). Do NOT recompute.

# ── Power-Law Fit ────────────────────────────────────────────────────────────
SLOPE_TARGET           = 0.0           # updated from 1/3 — density contrast prediction
FIT_METHOD             = "ODR"         # orthogonal distance regression (primary)
FIT_METHOD_SECONDARY   = "OLS"         # ordinary least squares (reported alongside)

# ── Provisional Thresholds (to be finalised by Millennium mock) ──────────────
# Source: Grok Checkpoint 2, 13 March 2026
# These values are PROVISIONAL. Final values set in finalise_thresholds()
# after mock calibration completes. Do not use for real-data decisions.
PROVISIONAL_R_OUTCOME_A     = 0.55     # Pearson r threshold for strong scaling
PROVISIONAL_R_OUTCOME_B_LO  = 0.35
PROVISIONAL_R_OUTCOME_B_HI  = 0.55
PROVISIONAL_CHI2_MAX        = 1.5
PROVISIONAL_SLOPE_TOL_STAT  = 0.08    # |α − 1/3| < 0.08 practical band
MIN_N_VOIDS                 = 800      # minimum sample after cuts

# Final thresholds locked after mock calibration (17 March 2026):
FINAL_R_OUTCOME_A        = -0.0732   # r must exceed this (5th pct null)
FINAL_R_OUTCOME_B_LO     = -0.0732   # 5th pct null
FINAL_R_OUTCOME_B_HI     =  0.0729   # 95th pct null
FINAL_MOCK_MEDIAN_R      =  0.0002
FINAL_MOCK_5PCT_R        = -0.0732
FINAL_MOCK_95PCT_R       =  0.0729
FINAL_MOCK_99PCT_R       =  0.0999
MOCK_CALIBRATION_DATE    = "2026-03-17"
MOCK_N_VOIDS             = 474
MOCK_SLOPE_ALPHA         = -0.0600
MOCK_SLOPE_SIGMA         =  0.0131
MOCK_SOURCE              = "Rodriguez-Medrano et al. 2024 TNG300 Popcorn voids"

# ── Null Test ────────────────────────────────────────────────────────────────
NULL_SHUFFLE_N         = 1000
NULL_BEAT_PERCENTILE   = 95.0          # real fit must beat 95% of shuffles
RANDOM_SEED            = 20260315      # date-stamped seed, do not change

# ── Simulation Comparison ────────────────────────────────────────────────────
SIM_NAME     = "IllustrisTNG TNG300-1 (public halo catalogues + void finder)"
SIM_URL      = "https://www.tng-project.org/data/downloads/TNG300-1/"
SIM_SNAPSHOT = 67          # z ≈ 0.5, centre of BOSS redshift window 0.2–0.7
SIM_VERSION  = "v1.1 — updated 15 March 2026, Millennium replaced by TNG300-1"
SIM_VOID_CATALOGUE_REF = (
    "Rodriguez-Medrano et al. (2024) MNRAS — TNG300 void catalogue; "
    "fallback: VIDE/ZOBOV applied to TNG300-1 snapshot 67 FoF groups"
)

# ── Falsification Statement (to be finalised at Checkpoint 3) ───────────────
# Placeholder — exact wording locked before real-data run
FALSIFICATION_STATEMENT = (
    "VOID-DOMAIN v1 falsification statement (locked 17 March 2026): "
    "If the fitted slope α deviates from 0 at >95% CL "
    "OR the Pearson r < -0.0732 (mock 5th percentile null), "
    "Outcome C is declared. This is a genuine FAIL for the "
    "gravitational-influence-domain fixed-contrast hypothesis. "
    "It does not falsify the broader TSM2 framework."
)

# ── Robustness Annulus Variants ──────────────────────────────────────────────
ROBUSTNESS_ANNULI = [
    (1.5, 3.5),    # inner sensitivity check
    (2.5, 5.0),    # outer sensitivity check
]

# ── Minimum Clusters per Void ────────────────────────────────────────────────
MIN_CLUSTERS_PER_VOID  = 3             # voids with fewer clusters excluded
