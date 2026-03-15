import os
import pickle
import urllib.request
import urllib.error
import numpy as np
import pandas as pd
from scipy import odr
from scipy.stats import pearsonr

from void_domain.manifest import (
    ANNULUS_INNER_FACTOR,
    ANNULUS_OUTER_FACTOR,
    MIN_CLUSTERS_PER_VOID,
    NULL_SHUFFLE_N,
    RANDOM_SEED,
    SIM_NAME,
    SIM_URL,
    SIM_SNAPSHOT,
    SIM_VERSION,
    MASS_UNIT,
)

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
CHECKPOINTS_DIR = os.path.join(os.path.dirname(__file__), "checkpoints")

GROUPS_PATH = os.path.join(DATA_DIR, "tng300_groups_snap67.hdf5")
GROUPS_CHECKPOINT = os.path.join(CHECKPOINTS_DIR, "tng300_groups_downloaded.pkl")
VOIDS_PATH_FITS = os.path.join(DATA_DIR, "tng300_voids.fits")
VOIDS_PATH_CSV = os.path.join(DATA_DIR, "tng300_voids.csv")
MOCK_TABLE_PATH = os.path.join(CHECKPOINTS_DIR, "mock_void_table.pkl")
FIT_RESULTS_PATH = os.path.join(CHECKPOINTS_DIR, "mock_fit_results.pkl")


def download_tng300_data():
    os.makedirs(DATA_DIR, exist_ok=True)
    os.makedirs(CHECKPOINTS_DIR, exist_ok=True)

    groups_ok = False
    groups_url = f"{SIM_URL}groups_067/"
    print(f"[1/2] Attempting TNG300-1 group catalogue download...")
    print(f"      URL: {groups_url}")
    try:
        req = urllib.request.Request(groups_url, method="HEAD")
        req.add_header("User-Agent", "VOID-DOMAIN/1.0")
        urllib.request.urlopen(req, timeout=15)
        print("      HEAD request succeeded — but full HDF5 download requires API key.")
        print("      Proceeding to fallback: synthetic mock will be used.")
    except (urllib.error.URLError, urllib.error.HTTPError, OSError) as e:
        print(f"      Download failed: {e}")
        print("      Proceeding to fallback.")

    voids_ok = False
    primary_url = "https://catalogs.iate.conicet.unc.edu.ar"
    print(f"\n[2/2] Attempting TNG300 void catalogue download...")
    print(f"      Primary URL: {primary_url}")
    try:
        req = urllib.request.Request(primary_url, method="HEAD")
        req.add_header("User-Agent", "VOID-DOMAIN/1.0")
        urllib.request.urlopen(req, timeout=15)
        print("      Catalogue server reachable — but specific TNG300 void file path unknown.")
        print("      VOID CATALOGUE NOT FOUND — switching to synthetic mock")
    except (urllib.error.URLError, urllib.error.HTTPError, OSError) as e:
        print(f"      Download failed: {e}")
        print("      VOID CATALOGUE NOT FOUND — switching to synthetic mock")

    return groups_ok, voids_ok


def _build_synthetic_mock():
    print("\nSYNTHETIC MOCK ACTIVE — real TNG300 data unavailable")
    rng = np.random.default_rng(RANDOM_SEED)

    n_voids = 1200

    log_r = rng.normal(loc=np.log10(35.0), scale=0.25, size=n_voids)
    r_void = 10.0 ** log_r

    log_m = rng.normal(loc=np.log10(5.0), scale=0.35, size=n_voids)
    m_surr = 10.0 ** log_m
    m_surr += rng.poisson(lam=2.0, size=n_voids) * 0.5

    n_groups = rng.poisson(lam=8, size=n_voids)
    z_void = rng.uniform(0.2, 0.7, size=n_voids)

    df = pd.DataFrame({
        "void_id": np.arange(n_voids),
        "R_void_hmpc": r_void,
        "M_surrounding_1e14Msun": m_surr,
        "z_void": z_void,
        "N_groups_in_annulus": n_groups,
    })

    df = df[df["N_groups_in_annulus"] >= MIN_CLUSTERS_PER_VOID].reset_index(drop=True)

    return df


def build_mock_void_table():
    if os.path.exists(GROUPS_PATH) and (os.path.exists(VOIDS_PATH_FITS) or os.path.exists(VOIDS_PATH_CSV)):
        print("Real TNG300 data files found — but parser not yet implemented.")
        print("Falling back to synthetic mock.")

    df = _build_synthetic_mock()

    os.makedirs(CHECKPOINTS_DIR, exist_ok=True)
    df.to_pickle(MOCK_TABLE_PATH)
    print(f"Mock void table built: N = {len(df)} voids after cuts")
    return df


def _odr_fit(x, y):
    def linear_func(B, x):
        return B[0] * x + B[1]

    model = odr.Model(linear_func)
    data = odr.RealData(x, y)
    initial = odr.ODR(data, model, beta0=[0.3, 0.0])
    output = initial.run()
    alpha = output.beta[0]
    sigma_alpha = output.sd_beta[0]
    return alpha, sigma_alpha


def run_mock_fit(df):
    log_r = np.log10(df["R_void_hmpc"].values)
    log_m = np.log10(df["M_surrounding_1e14Msun"].values)

    alpha, sigma_alpha = _odr_fit(log_m, log_r)
    r_val, _ = pearsonr(log_m, log_r)
    n = len(df)
    residuals = log_r - (alpha * log_m + 0.0)
    chi2_dof = np.sum(residuals ** 2) / (n - 2)

    print(f"\n--- Mock real fit ---")
    print(f"  α = {alpha:.4f} ± {sigma_alpha:.4f}")
    print(f"  Pearson r = {r_val:.4f}")
    print(f"  Reduced χ²/dof = {chi2_dof:.4f}")

    rng = np.random.default_rng(RANDOM_SEED)
    null_r_dist = np.empty(NULL_SHUFFLE_N)
    null_alpha_dist = np.empty(NULL_SHUFFLE_N)

    print(f"\nRunning {NULL_SHUFFLE_N} null shuffles...")
    for i in range(NULL_SHUFFLE_N):
        shuffled_m = rng.permutation(log_m)
        a_shuf, _ = _odr_fit(shuffled_m, log_r)
        r_shuf, _ = pearsonr(shuffled_m, log_r)
        null_r_dist[i] = r_shuf
        null_alpha_dist[i] = a_shuf

    null_r_median = np.median(null_r_dist)
    null_r_95pct = np.percentile(null_r_dist, 95)
    null_r_99pct = np.percentile(null_r_dist, 99)

    results = {
        "alpha": alpha,
        "sigma_alpha": sigma_alpha,
        "mock_real_r": r_val,
        "chi2_dof": chi2_dof,
        "null_r_dist": null_r_dist,
        "null_alpha_dist": null_alpha_dist,
        "null_r_median": null_r_median,
        "null_r_95pct": null_r_95pct,
        "null_r_99pct": null_r_99pct,
        "n_voids": n,
    }

    os.makedirs(CHECKPOINTS_DIR, exist_ok=True)
    with open(FIT_RESULTS_PATH, "wb") as f:
        pickle.dump(results, f)

    return results


def finalise_thresholds(results):
    null_r_95pct = results["null_r_95pct"]

    final_r_outcome_a = null_r_95pct + 0.10
    final_r_outcome_b_lo = null_r_95pct
    final_r_outcome_b_hi = null_r_95pct + 0.10

    print(f"""
=== VOID-DOMAIN v1 Mock Calibration Report ===
Mock source:          {SIM_NAME}
N voids (mock):       {results['n_voids']}
Real mock fit:        α = {results['alpha']:.4f} ± {results['sigma_alpha']:.4f}, r = {results['mock_real_r']:.4f}
Null r (median):      {results['null_r_median']:.4f}
Null r (95th pct):    {results['null_r_95pct']:.4f}
Null r (99th pct):    {results['null_r_99pct']:.4f}
-------------------------------------------------
PROPOSED THRESHOLDS (pending Graham/Grok review):
  Outcome A (strong):  r > {final_r_outcome_a:.4f}
  Outcome B (weak):    {final_r_outcome_b_lo:.4f} < r ≤ {final_r_outcome_b_hi:.4f}
  Outcome C (FAIL):    r ≤ {final_r_outcome_b_lo:.4f}
================================================""")

    return {
        "FINAL_R_OUTCOME_A": final_r_outcome_a,
        "FINAL_R_OUTCOME_B_LO": final_r_outcome_b_lo,
        "FINAL_R_OUTCOME_B_HI": final_r_outcome_b_hi,
    }
