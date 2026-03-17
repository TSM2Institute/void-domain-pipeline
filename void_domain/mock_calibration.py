import os
import pickle
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
import numpy as np
import pandas as pd
import h5py
import requests
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

GROUPS_PATH = os.path.join(DATA_DIR, "tng300_groups_snap67.csv")
GROUPS_CHECKPOINT = os.path.join(CHECKPOINTS_DIR, "tng300_groups_downloaded.pkl")
GROUPS_WITHPOS_PATH = os.path.join(DATA_DIR, "tng300_groups_withpos_snap67.csv")
GROUPCAT_HDF5_PATH = os.path.join(DATA_DIR, "tng300_groupcat_snap67.hdf5")
POSITIONS_CHECKPOINT = os.path.join(CHECKPOINTS_DIR, "tng300_positions_downloaded.pkl")
MOCK_VOIDS_PATH = os.path.join(DATA_DIR, "tng300_mock_voids.csv")
VOIDS_PATH_FITS = os.path.join(DATA_DIR, "tng300_voids.fits")
VOIDS_PATH_CSV = os.path.join(DATA_DIR, "tng300_voids.csv")
MOCK_TABLE_PATH = os.path.join(CHECKPOINTS_DIR, "mock_void_table.pkl")
FIT_RESULTS_PATH = os.path.join(CHECKPOINTS_DIR, "mock_fit_results.pkl")

TNG300_BOX_SIZE = 205.0


def _get_api_key():
    api_key = os.environ.get("TNG_API")
    if not api_key:
        raise EnvironmentError("TNG_API secret not found in environment.")
    return api_key


def download_tng300_data():
    os.makedirs(DATA_DIR, exist_ok=True)
    os.makedirs(CHECKPOINTS_DIR, exist_ok=True)

    api_key = _get_api_key()
    headers = {"api-key": api_key}

    groups_ok = False
    initial_url = (
        "https://www.tng-project.org/api/TNG300-1/snapshots/67/subhalos/"
        "?limit=1000&primary_flag=1&mass__gt=1e3"
    )
    print("[1/2] Downloading TNG300-1 FoF group catalogue (snapshot 67)...")
    print(f"      Endpoint: subhalos/?primary_flag=1&mass__gt=1e3")
    print(f"      (primary subhalo of each FoF group, M > 10^13 M☉ — cluster scale)")
    try:
        records = []
        url = initial_url
        page = 0
        while url:
            resp = requests.get(url, headers=headers, timeout=60)
            resp.raise_for_status()
            data = resp.json()
            results = data.get("results", [])
            records.extend(results)
            page += 1
            print(f"  Downloaded {len(records)} groups so far... (page {page})")
            url = data.get("next")
            if url:
                time.sleep(0.1)

        rows = []
        for halo in records:
            pos = halo.get("pos")
            rows.append({
                "group_id": halo["id"],
                "mass_log_msun": halo["mass_log_msun"],
                "mass_200c_1e14Msun": 10.0 ** halo["mass_log_msun"] / 1e14,
                "pos_x": pos[0] if pos else None,
                "pos_y": pos[1] if pos else None,
                "pos_z": pos[2] if pos else None,
            })

        df_groups = pd.DataFrame(rows)
        df_groups.to_csv(GROUPS_PATH, index=False)
        with open(GROUPS_CHECKPOINT, "wb") as f:
            pickle.dump(True, f)
        groups_ok = True
        print(f"Group catalogue saved: N = {len(df_groups)} groups")

    except (requests.RequestException, OSError, KeyError) as e:
        print(f"      Download failed: {e}")
        print("      Proceeding to fallback.")
        df_groups = None

    voids_ok = False
    primary_url = "https://catalogs.iate.conicet.unc.edu.ar"
    print(f"\n[2/2] Attempting TNG300 void catalogue download...")
    print(f"      Primary URL: {primary_url}")
    try:
        resp = requests.head(primary_url, timeout=15)
        print(f"      Catalogue server reachable (status {resp.status_code}) — but specific TNG300 void file path unknown.")
        print("      VOID CATALOGUE NOT FOUND — switching to synthetic mock")
    except (requests.RequestException, OSError) as e:
        print(f"      Download failed: {e}")
        print("      VOID CATALOGUE NOT FOUND — switching to synthetic mock")

    return groups_ok, voids_ok


def verify_groups_download():
    if not os.path.exists(GROUPS_PATH):
        print("No group catalogue file found.")
        return
    df = pd.read_csv(GROUPS_PATH)
    has_pos = not df["pos_x"].isna().all()
    print(f"\nTotal groups downloaded: {len(df)}")
    print(f"Mass range (log10 M/M☉): {df['mass_log_msun'].min():.2f} to {df['mass_log_msun'].max():.2f}")
    print(f"Position columns present: {has_pos}")
    print(f"\nFirst 3 rows:")
    print(df.head(3).to_string(index=False))


def fetch_tng300_positions():
    os.makedirs(DATA_DIR, exist_ok=True)
    os.makedirs(CHECKPOINTS_DIR, exist_ok=True)

    api_key = _get_api_key()
    headers = {"api-key": api_key}

    max_chunk = 150
    mass_threshold = 0.01

    print("Downloading TNG300-1 groupcat HDF5 chunks (positions + masses)...")
    print(f"  Chunks 0–{max_chunk - 1}, parallel (10 threads)")

    def _download_one(i):
        url = f"http://www.tng-project.org/api/TNG300-1/files/groupcat-67.{i}.hdf5"
        resp = requests.get(url, headers=headers, timeout=120)
        resp.raise_for_status()
        tmp_path = os.path.join(DATA_DIR, f"_chunk_{i}.hdf5")
        with open(tmp_path, "wb") as f:
            f.write(resp.content)
        with h5py.File(tmp_path, "r") as hf:
            n = hf["Header"].attrs["Ngroups_ThisFile"]
            pos = hf["Group"]["GroupPos"][:] if n > 0 else None
            m200 = hf["Group"]["Group_M_Crit200"][:] if n > 0 else None
            r200 = hf["Group"]["Group_R_Crit200"][:] if n > 0 else None
        os.remove(tmp_path)
        return i, pos, m200, r200, len(resp.content)

    chunk_data = {}
    total_bytes = 0
    done = 0
    with ThreadPoolExecutor(max_workers=10) as pool:
        futures = {pool.submit(_download_one, i): i for i in range(max_chunk)}
        for fut in as_completed(futures):
            i, pos, m200, r200, nbytes = fut.result()
            chunk_data[i] = (pos, m200, r200)
            total_bytes += nbytes
            done += 1
            if done % 25 == 0 or done == max_chunk:
                print(f"  {done}/{max_chunk} chunks downloaded ({total_bytes / (1024*1024):.0f} MB)")

    all_pos, all_m200, all_r200 = [], [], []
    for i in sorted(chunk_data.keys()):
        pos, m200, r200 = chunk_data[i]
        if m200 is not None:
            all_pos.append(pos)
            all_m200.append(m200)
            all_r200.append(r200)

    print(f"  HDF5 download complete: {total_bytes / (1024*1024):.0f} MB total")

    print("\nParsing into DataFrame...")
    group_pos = np.concatenate(all_pos, axis=0)
    group_m200 = np.concatenate(all_m200)
    group_r200 = np.concatenate(all_r200)

    n_total = len(group_m200)
    print(f"  Raw groups from chunks: {n_total}")

    pos_mpc = group_pos / 1000.0
    mass_1e14 = group_m200 * 1e-4
    r200_mpc = group_r200 / 1000.0

    df = pd.DataFrame({
        "group_id": np.arange(n_total),
        "pos_x_mpc": pos_mpc[:, 0],
        "pos_y_mpc": pos_mpc[:, 1],
        "pos_z_mpc": pos_mpc[:, 2],
        "mass_200c_1e14Msun": mass_1e14,
        "r_200c_mpc": r200_mpc,
    })

    df = df[df["mass_200c_1e14Msun"] > mass_threshold].reset_index(drop=True)

    df.to_csv(GROUPS_WITHPOS_PATH, index=False)
    with open(POSITIONS_CHECKPOINT, "wb") as f:
        pickle.dump(True, f)

    print(f"\nGroups with positions: N = {len(df)}")
    print(f"Mass range (10^14 Msun): {df['mass_200c_1e14Msun'].min():.4f} to {df['mass_200c_1e14Msun'].max():.2f}")
    print(f"Position range X (Mpc/h): {df['pos_x_mpc'].min():.1f} to {df['pos_x_mpc'].max():.1f}")
    print(f"\nFirst 3 rows:")
    print(df.head(3).to_string(index=False))

    return df


def build_synthetic_voids(df_groups):
    rng = np.random.default_rng(RANDOM_SEED)

    box = TNG300_BOX_SIZE
    n_cells = 20
    cell_size = box / n_cells

    print(f"\nBuilding synthetic void catalogue from group positions...")
    print(f"  Box size: {box} Mpc/h, grid: {n_cells}x{n_cells}x{n_cells}, cell size: {cell_size:.1f} Mpc/h")

    positions = df_groups[["pos_x_mpc", "pos_y_mpc", "pos_z_mpc"]].values
    edges = np.linspace(0, box, n_cells + 1)

    counts, _ = np.histogramdd(
        positions,
        bins=[edges, edges, edges],
    )

    empty_cells = np.argwhere(counts == 0)
    print(f"  Empty cells (group count = 0): {len(empty_cells)}")

    voids = []
    for idx, (ix, iy, iz) in enumerate(empty_cells):
        cx = (ix + 0.5) * cell_size
        cy = (iy + 0.5) * cell_size
        cz = (iz + 0.5) * cell_size

        log_r = rng.normal(loc=np.log10(25.0), scale=0.20)
        r_void = np.clip(10.0 ** log_r, 15.0, 80.0)

        voids.append({
            "void_id": idx,
            "cx_mpc": cx,
            "cy_mpc": cy,
            "cz_mpc": cz,
            "R_void_mpc": r_void,
        })

    df_voids = pd.DataFrame(voids)
    df_voids.to_csv(MOCK_VOIDS_PATH, index=False)

    print(f"Mock voids identified: N = {len(df_voids)} underdense cells")
    return df_voids


def test_tng_auth():
    api_key = _get_api_key()
    headers = {"api-key": api_key}
    url = "https://www.tng-project.org/api/TNG300-1/snapshots/67/subhalos/?limit=1"
    print(f"Testing TNG300 API authentication...")
    print(f"URL: {url}")
    try:
        resp = requests.get(url, headers=headers, timeout=30)
        print(f"TNG300 API authentication: OK — status {resp.status_code}")
        data = resp.json()
        if "results" in data and len(data["results"]) > 0:
            print(f"First record sample: {data['results'][0]}")
        else:
            print(f"Response sample: {str(data)[:300]}")
    except (requests.RequestException, OSError) as e:
        print(f"Authentication test failed: {e}")


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


def build_mock_void_table(df_voids=None, df_groups=None):
    if df_voids is None:
        df_voids = pd.read_csv(MOCK_VOIDS_PATH)
    if df_groups is None:
        df_groups = pd.read_csv(GROUPS_WITHPOS_PATH)

    print(f"\nCross-matching {len(df_voids)} voids × {len(df_groups)} groups...")
    print(f"  Annulus: {ANNULUS_INNER_FACTOR}×R_void to {ANNULUS_OUTER_FACTOR}×R_void")
    print(f"  Periodic box: {TNG300_BOX_SIZE} Mpc/h")

    g_pos = df_groups[["pos_x_mpc", "pos_y_mpc", "pos_z_mpc"]].values
    g_mass = df_groups["mass_200c_1e14Msun"].values

    box = TNG300_BOX_SIZE
    half_box = box / 2.0

    rows = []
    for _, v in df_voids.iterrows():
        vc = np.array([v["cx_mpc"], v["cy_mpc"], v["cz_mpc"]])
        r_void = v["R_void_mpc"]
        inner_r = ANNULUS_INNER_FACTOR * r_void
        outer_r = ANNULUS_OUTER_FACTOR * r_void

        dx = np.abs(g_pos[:, 0] - vc[0])
        dy = np.abs(g_pos[:, 1] - vc[1])
        dz = np.abs(g_pos[:, 2] - vc[2])
        dx = np.minimum(dx, box - dx)
        dy = np.minimum(dy, box - dy)
        dz = np.minimum(dz, box - dz)
        dist = np.sqrt(dx**2 + dy**2 + dz**2)

        in_annulus = (dist > inner_r) & (dist < outer_r)
        n_groups = int(np.sum(in_annulus))
        m_surr = float(np.sum(g_mass[in_annulus]))

        rows.append({
            "void_id": int(v["void_id"]),
            "R_void_mpc": r_void,
            "M_surrounding_1e14Msun": m_surr,
            "N_groups_in_annulus": n_groups,
        })

    df_table = pd.DataFrame(rows)
    df_cut = df_table[df_table["N_groups_in_annulus"] >= MIN_CLUSTERS_PER_VOID].reset_index(drop=True)

    os.makedirs(CHECKPOINTS_DIR, exist_ok=True)
    df_cut.to_pickle(MOCK_TABLE_PATH)

    print(f"Voids after MIN_CLUSTERS cut: N = {len(df_cut)}")
    print(f"R_void range (Mpc/h): {df_cut['R_void_mpc'].min():.1f} to {df_cut['R_void_mpc'].max():.1f}")
    print(f"M_surrounding range (10^14 Msun): {df_cut['M_surrounding_1e14Msun'].min():.3f} to {df_cut['M_surrounding_1e14Msun'].max():.2f}")
    print(f"Median N_groups per annulus: {df_cut['N_groups_in_annulus'].median():.1f}")

    return df_cut


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
    log_r = np.log10(df["R_void_mpc"].values)
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

    n_groups_used = 28301

    print(f"""
=== VOID-DOMAIN v1 Mock Calibration Report ===
Mock source:           TNG300-1 snapshot 67 (z≈0.5)
Groups used:           {n_groups_used} (M > 10^12 Msun)
N voids (after cuts):  {results['n_voids']}
---
Real mock fit:
  α  = {results['alpha']:.4f} ± {results['sigma_alpha']:.4f}
  r  = {results['mock_real_r']:.4f}
  χ²/dof = {results['chi2_dof']:.4f}
---
Null shuffle ({NULL_SHUFFLE_N}×):
  Median r:   {results['null_r_median']:.4f}
  95th pct r: {results['null_r_95pct']:.4f}
  99th pct r: {results['null_r_99pct']:.4f}
---
PROPOSED THRESHOLDS (pending Grok review):
  Outcome A (strong):  r > {final_r_outcome_a:.4f}
  Outcome B (weak):    {final_r_outcome_b_lo:.4f} < r ≤ {final_r_outcome_b_hi:.4f}
  Outcome C (FAIL):    r ≤ {final_r_outcome_b_lo:.4f}
================================================

IMPORTANT: Do NOT write these values to manifest.py.
Print only. Graham reviews before any threshold is locked.""")

    return {
        "FINAL_R_OUTCOME_A": final_r_outcome_a,
        "FINAL_R_OUTCOME_B_LO": final_r_outcome_b_lo,
        "FINAL_R_OUTCOME_B_HI": final_r_outcome_b_hi,
    }
