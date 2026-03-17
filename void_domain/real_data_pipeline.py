import os
import pickle
import numpy as np
import pandas as pd
from astropy.cosmology import FlatLambdaCDM

from scipy import stats as sp_stats
from scipy.odr import ODR, Model, RealData

from void_domain.manifest import (
    ANNULUS_INNER_FACTOR,
    ANNULUS_OUTER_FACTOR,
    MIN_CLUSTERS_PER_VOID,
    VOID_Z_MIN,
    VOID_Z_MAX,
    FINAL_MOCK_5PCT_R,
    FALSIFICATION_STATEMENT,
    RANDOM_SEED,
)

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
CHECKPOINTS_DIR = os.path.join(os.path.dirname(__file__), "checkpoints")
OUTPUTS_DIR = os.path.join(os.path.dirname(__file__), "outputs")

MAO_VOIDS_PATH = os.path.join(DATA_DIR, "mao2017_boss_voids.dat")
REDMAPPER_PATH = os.path.join(DATA_DIR, "redmapper_dr8_vizier.dat")
REAL_TABLE_PKL = os.path.join(CHECKPOINTS_DIR, "real_void_table.pkl")
REAL_TABLE_CSV = os.path.join(OUTPUTS_DIR, "real_void_table.csv")

COSMO = FlatLambdaCDM(H0=70, Om0=0.3)

BOSS_SKY_FRACTION = 0.2424

Z_BINS = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]

LAMBDA_MIN = 20.0
CLUSTER_Z_MIN = 0.1
CLUSTER_Z_MAX = 0.6


def load_voids():
    col_names = [
        "Sample", "ID", "RAdeg", "DEdeg", "z", "NGal",
        "V", "Reff", "nmin", "delmin", "r", "Prob", "Dbound",
    ]
    df = pd.read_fwf(
        MAO_VOIDS_PATH,
        colspecs=[
            (0, 11), (12, 17), (18, 25), (26, 32), (33, 38),
            (39, 45), (46, 55), (56, 63), (64, 73), (74, 80),
            (81, 86), (87, 96), (97, 104),
        ],
        names=col_names,
    )
    df = df[(df["z"] > VOID_Z_MIN) & (df["z"] < VOID_Z_MAX)].reset_index(drop=True)
    print(f"Voids after redshift cut: N = {len(df)}")
    return df


def load_clusters():
    df = pd.read_csv(REDMAPPER_PATH, sep="\t", comment="#", low_memory=False)
    df = df[~df.iloc[:, 0].astype(str).str.match(r"^\s*$|^\s*-+")].reset_index(drop=True)
    df = df.iloc[1:].reset_index(drop=True)

    df["zlambda"] = pd.to_numeric(df["zlambda"], errors="coerce")
    df["lambda"] = pd.to_numeric(df["lambda"], errors="coerce")
    df["RAJ2000"] = pd.to_numeric(df["RAJ2000"], errors="coerce")
    df["DEJ2000"] = pd.to_numeric(df["DEJ2000"], errors="coerce")

    df = df.dropna(subset=["zlambda", "lambda", "RAJ2000", "DEJ2000"])
    df = df[(df["zlambda"] > CLUSTER_Z_MIN) & (df["zlambda"] < CLUSTER_Z_MAX)].reset_index(drop=True)
    df = df[df["lambda"] >= LAMBDA_MIN].reset_index(drop=True)
    print(f"Clusters after cuts: N = {len(df)}")
    return df


def radec_to_cartesian(ra_deg, dec_deg, z):
    d_c = COSMO.comoving_distance(z).value
    ra_rad = np.deg2rad(ra_deg)
    dec_rad = np.deg2rad(dec_deg)
    x = d_c * np.cos(dec_rad) * np.cos(ra_rad)
    y = d_c * np.cos(dec_rad) * np.sin(ra_rad)
    z_cart = d_c * np.sin(dec_rad)
    return x, y, z_cart


def lambda_to_m200c(lam):
    log_m = 14.344 + 1.33 * np.log10(lam / 40.0)
    return 10.0 ** log_m / 1e14


def compute_redshift_binned_density(df_clusters, mass_1e14):
    print("\nRedshift-binned mean densities (10^14 Msun/Mpc^3):")
    rho_bins = {}
    for i in range(len(Z_BINS) - 1):
        z_lo, z_hi = Z_BINS[i], Z_BINS[i + 1]
        z_mid = (z_lo + z_hi) / 2.0
        mask = (df_clusters["zlambda"].values >= z_lo) & (df_clusters["zlambda"].values < z_hi)
        total_mass = mass_1e14[mask].sum()
        d_lo = COSMO.comoving_distance(z_lo).value
        d_hi = COSMO.comoving_distance(z_hi).value
        v_shell = (4.0 / 3.0) * np.pi * (d_hi ** 3 - d_lo ** 3)
        v_eff = BOSS_SKY_FRACTION * v_shell
        rho = total_mass / v_eff if v_eff > 0 else 0.0
        rho_bins[z_mid] = rho
        n_cl = int(mask.sum())
        print(f"  z={z_mid:.2f}: rho={rho:.6f}, N_clusters={n_cl}, V_eff={v_eff:.0f} Mpc^3")
    return rho_bins


def get_rho_for_z(z_void, rho_bins):
    bin_mids = np.array(sorted(rho_bins.keys()))
    idx = np.argmin(np.abs(bin_mids - z_void))
    return rho_bins[bin_mids[idx]]


def run():
    os.makedirs(CHECKPOINTS_DIR, exist_ok=True)
    os.makedirs(OUTPUTS_DIR, exist_ok=True)

    print("=" * 60)
    print("VOID-DOMAIN v1 — Real Data Cross-Match Pipeline")
    print("=" * 60)

    print("\n--- STEP 1: Load and filter catalogues ---")
    df_voids = load_voids()
    df_clusters = load_clusters()

    print("\n--- STEP 2: Convert to comoving Cartesian ---")
    vx, vy, vz_cart = radec_to_cartesian(
        df_voids["RAdeg"].values,
        df_voids["DEdeg"].values,
        df_voids["z"].values,
    )
    print(f"Void positions computed: d_c range {np.sqrt(vx**2+vy**2+vz_cart**2).min():.0f} to {np.sqrt(vx**2+vy**2+vz_cart**2).max():.0f} Mpc")

    cx, cy, cz_cart = radec_to_cartesian(
        df_clusters["RAJ2000"].values,
        df_clusters["DEJ2000"].values,
        df_clusters["zlambda"].values,
    )
    print(f"Cluster positions computed: d_c range {np.sqrt(cx**2+cy**2+cz_cart**2).min():.0f} to {np.sqrt(cx**2+cy**2+cz_cart**2).max():.0f} Mpc")

    print("\n--- STEP 3: Convert lambda to M200c ---")
    lam_vals = df_clusters["lambda"].values
    mass_1e14 = lambda_to_m200c(lam_vals)
    print(f"Mass range (10^14 Msun/h): {mass_1e14.min():.3f} to {mass_1e14.max():.2f}")

    print("\n--- STEP 4: Redshift-binned mean density ---")
    rho_bins = compute_redshift_binned_density(df_clusters, mass_1e14)

    print("\n--- STEP 5: Cross-match (annulus 2-4×R_void) ---")
    c_pos = np.column_stack([cx, cy, cz_cart])

    rows = []
    n_voids = len(df_voids)
    for i in range(n_voids):
        r_void = df_voids.iloc[i]["Reff"]
        inner_r = ANNULUS_INNER_FACTOR * r_void
        outer_r = ANNULUS_OUTER_FACTOR * r_void

        dx = c_pos[:, 0] - vx[i]
        dy = c_pos[:, 1] - vy[i]
        dz = c_pos[:, 2] - vz_cart[i]
        dist = np.sqrt(dx ** 2 + dy ** 2 + dz ** 2)

        in_ann = (dist > inner_r) & (dist < outer_r)
        n_in = int(in_ann.sum())
        m_ann = float(mass_1e14[in_ann].sum())

        v_ann = (4.0 / 3.0) * np.pi * 56.0 * r_void ** 3
        rho_ann = m_ann / v_ann if v_ann > 0 else 0.0
        rho_mean = get_rho_for_z(df_voids.iloc[i]["z"], rho_bins)
        delta = (rho_ann - rho_mean) / rho_mean if rho_mean > 0 else 0.0

        rows.append({
            "void_id": int(df_voids.iloc[i]["ID"]),
            "sample": df_voids.iloc[i]["Sample"],
            "RA": df_voids.iloc[i]["RAdeg"],
            "Dec": df_voids.iloc[i]["DEdeg"],
            "z_void": df_voids.iloc[i]["z"],
            "R_void_mpc": r_void,
            "M_annulus_1e14": m_ann,
            "N_in_annulus": n_in,
            "V_annulus_mpc3": v_ann,
            "rho_annulus": rho_ann,
            "rho_mean_bin": rho_mean,
            "delta": delta,
        })

        if (i + 1) % 200 == 0 or i == n_voids - 1:
            print(f"  {i + 1}/{n_voids} voids processed")

    df_result = pd.DataFrame(rows)

    n_before = len(df_result)
    df_cut = df_result[df_result["N_in_annulus"] >= MIN_CLUSTERS_PER_VOID].reset_index(drop=True)
    n_after = len(df_cut)

    df_cut.to_pickle(REAL_TABLE_PKL)
    df_cut.to_csv(REAL_TABLE_CSV, index=False)

    n_high_z = len(df_voids[df_voids["z"] > 0.55])

    print(f"\n--- STEP 6: Cross-Match Report ---")
    print(f"""
=== VOID-DOMAIN v1 Real Data Cross-Match Report ===
Void catalogue:    Mao et al. (2017) BOSS DR12
Cluster catalogue: redMaPPer DR8 (Rykoff et al. 2014)
Redshift note:     Clusters only to z=0.55; voids to z=0.70
---
Voids after z cut:      {n_before}
Clusters after cuts:    {len(df_clusters)}
Voids after MIN_CLUSTERS cut: {n_after}
Voids lost (z>0.55):    {n_high_z} (estimated)
---
R_void range (Mpc/h):      {df_cut['R_void_mpc'].min():.1f} to {df_cut['R_void_mpc'].max():.1f}
delta range:               {df_cut['delta'].min():.4f} to {df_cut['delta'].max():.4f}
Median delta:              {df_cut['delta'].median():.4f}
Negative delta voids:      {(df_cut['delta'] < 0).sum()} / {n_after}
Median N_clusters/annulus: {df_cut['N_in_annulus'].median():.1f}
---
RAW DATA ONLY — fit and outcome pending Grok CP3 review
====================================================""")


REAL_CUT_PKL = os.path.join(CHECKPOINTS_DIR, "real_void_table_cut.pkl")
REAL_FIT_PKL = os.path.join(CHECKPOINTS_DIR, "real_fit_results.pkl")


def run_real_fit():
    df = pd.read_pickle(REAL_CUT_PKL)
    n_all = len(df)

    df_pos = df[df["delta"] > 0].reset_index(drop=True)
    n_pos = len(df_pos)
    n_neg = n_all - n_pos

    print(f"Voids with delta > 0 (log-log fit): N = {n_pos}")
    print(f"Voids with delta <= 0 (excluded from log-log): N = {n_neg}")

    print("\n--- PART A: Linear fit (all voids) ---")
    delta_all = df["delta"].values
    r_all = df["R_void_mpc"].values
    slope, intercept, r_linear, p_linear, se = sp_stats.linregress(delta_all, r_all)
    print("Linear fit (all voids):")
    print(f"  slope = {slope:.4f} ± {se:.4f}")
    print(f"  Pearson r = {r_linear:.4f}")
    print(f"  p-value = {p_linear:.4e}")

    print("\n--- PART B: Log-log fit (delta > 0 voids) ---")
    x = np.log10(df_pos["delta"].values)
    y = np.log10(df_pos["R_void_mpc"].values)

    def linear_func(B, x):
        return B[0] * x + B[1]

    model = Model(linear_func)
    data = RealData(x, y)
    odr_obj = ODR(data, model, beta0=[0.0, 1.5])
    output = odr_obj.run()
    alpha = output.beta[0]
    beta_intercept = output.beta[1]
    sigma_alpha = output.sd_beta[0]
    r_loglog = np.corrcoef(x, y)[0, 1]
    residuals = y - (alpha * x + beta_intercept)
    chi2_dof = np.sum(residuals ** 2) / (len(x) - 2)
    print("Log-log fit (delta > 0 voids):")
    print(f"  α = {alpha:.4f} ± {sigma_alpha:.4f}")
    print(f"  Pearson r = {r_loglog:.4f}")
    print(f"  reduced χ²/dof = {chi2_dof:.4f}")
    print(f"  N used = {n_pos}")

    print("\n--- PART C: Bootstrap 95% CI on α ---")
    np.random.seed(RANDOM_SEED)
    alpha_boot = []
    for _ in range(1000):
        idx = np.random.choice(len(x), len(x), replace=True)
        xb = x[idx]
        yb = y[idx]
        data_b = RealData(xb, yb)
        odr_b = ODR(data_b, model, beta0=[0.0, 1.5])
        out_b = odr_b.run()
        alpha_boot.append(out_b.beta[0])
    alpha_boot = np.array(alpha_boot)
    ci_lo = np.percentile(alpha_boot, 2.5)
    ci_hi = np.percentile(alpha_boot, 97.5)
    print(f"Bootstrap 95% CI on α: [{ci_lo:.4f}, {ci_hi:.4f}]")
    alpha_zero_in_ci = (ci_lo <= 0.0 <= ci_hi)
    print(f"α = 0 inside 95% CI: {alpha_zero_in_ci}")

    print("\n--- PART D: 1000× label shuffle null test ---")
    np.random.seed(RANDOM_SEED)
    null_r = []
    delta_shuf = delta_all.copy()
    for _ in range(1000):
        np.random.shuffle(delta_shuf)
        _, _, r_s, _, _ = sp_stats.linregress(delta_shuf, r_all)
        null_r.append(r_s)
    null_r = np.array(null_r)
    p_null = np.mean(np.abs(null_r) >= np.abs(r_linear))
    print(f"Null shuffle p-value: {p_null:.4f}")
    print(f"Real |r| beats {100 * (1 - p_null):.1f}% of null shuffles")

    print("\n--- PART E: Decision Tree ---")
    print(f"Falsification statement:\n  {FALSIFICATION_STATEMENT}\n")

    r_above_floor = (r_linear > FINAL_MOCK_5PCT_R)

    if alpha_zero_in_ci and r_above_floor:
        outcome = "OUTCOME A — Strong support for fixed-contrast domains"
    elif not alpha_zero_in_ci and not r_above_floor:
        outcome = "OUTCOME C — FAIL (falsification statement fires)"
    elif not alpha_zero_in_ci or not r_above_floor:
        outcome = "OUTCOME C — FAIL (falsification statement fires)"
    else:
        outcome = "OUTCOME B — Weak/marginal"

    falsification_fires = not (alpha_zero_in_ci and r_above_floor)

    print("=== VOID-DOMAIN v1 FINAL RESULT ===")
    print(f"α = 0 inside 95% CI:  {alpha_zero_in_ci}")
    print(f"r > mock 5th pct ({FINAL_MOCK_5PCT_R}): {r_above_floor}")
    print(f"Falsification fires:  {falsification_fires}")
    print()
    print(f"*** {outcome} ***")
    print("====================================")

    results = {
        "slope_linear": slope,
        "se_linear": se,
        "r_linear": r_linear,
        "p_linear": p_linear,
        "alpha_loglog": alpha,
        "sigma_alpha": sigma_alpha,
        "r_loglog": r_loglog,
        "chi2_dof": chi2_dof,
        "n_pos": n_pos,
        "ci_lo": ci_lo,
        "ci_hi": ci_hi,
        "alpha_zero_in_ci": alpha_zero_in_ci,
        "null_r": null_r,
        "p_null": p_null,
        "r_above_floor": r_above_floor,
        "outcome": outcome,
    }
    with open(REAL_FIT_PKL, "wb") as f:
        pickle.dump(results, f)

    return results


if __name__ == "__main__":
    run()
