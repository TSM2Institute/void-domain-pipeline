import os
import numpy as np
import pandas as pd
import requests

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")

MAO_VOIDS_PATH = os.path.join(DATA_DIR, "mao2017_boss_voids.dat")
MAO_README_PATH = os.path.join(DATA_DIR, "mao2017_ReadMe.txt")
REDMAPPER_PATH = os.path.join(DATA_DIR, "redmapper_dr8_vizier.dat")
REDMAPPER_README_PATH = os.path.join(DATA_DIR, "redmapper_dr8_ReadMe.txt")


def download_file(urls, save_path, chunk_size=1024 * 1024, progress_interval_mb=10):
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    for url in urls:
        try:
            resp = requests.get(url, stream=True, timeout=120)
            resp.raise_for_status()
            total = 0
            with open(save_path, "wb") as f:
                for chunk in resp.iter_content(chunk_size=chunk_size):
                    f.write(chunk)
                    total += len(chunk)
                    mb = total / (1024 * 1024)
                    if progress_interval_mb and int(mb) > 0 and int(mb) % progress_interval_mb == 0:
                        prev = (total - len(chunk)) / (1024 * 1024)
                        if int(prev) % progress_interval_mb != 0 or int(prev) != int(mb):
                            print(f"  Downloaded {int(mb)} MB...")
            return total, url
        except (requests.RequestException, OSError) as e:
            print(f"  FAILED: {url} — {e}")
    return None, None


def download_vizier_tsv(vizier_table, save_path, timeout=120):
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    url = "https://vizier.cds.unistra.fr/viz-bin/asu-tsv"
    params = {
        "-source": vizier_table,
        "-out.max": "unlimited",
        "-out": "**",
    }
    try:
        resp = requests.get(url, params=params, timeout=timeout)
        resp.raise_for_status()
        with open(save_path, "wb") as f:
            f.write(resp.content)
        return len(resp.content)
    except (requests.RequestException, OSError) as e:
        print(f"  FAILED VizieR query for {vizier_table}: {e}")
        return None


def run():
    print("=" * 60)
    print("VOID-DOMAIN v1 — Real Data Acquisition")
    print("=" * 60)

    print("\n--- STEP 1: Mao et al. (2017) BOSS DR12 void catalogue ---")
    nbytes, _ = download_file(
        [
            "https://cdsarc.cds.unistra.fr/ftp/J/ApJ/835/161/table1.dat",
            "https://cdsarc.u-strasbg.fr/ftp/J/ApJ/835/161/table1.dat",
        ],
        MAO_VOIDS_PATH,
        progress_interval_mb=0,
    )
    if nbytes:
        print(f"Mao 2017 void catalogue downloaded: {nbytes // 1024} KB")
    else:
        print("ERROR: All Mao 2017 void catalogue URLs failed.")

    print("\n--- STEP 2: Mao 2017 ReadMe ---")
    nbytes_rm, _ = download_file(
        [
            "https://cdsarc.cds.unistra.fr/ftp/J/ApJ/835/161/ReadMe",
            "https://cdsarc.u-strasbg.fr/ftp/J/ApJ/835/161/ReadMe",
        ],
        MAO_README_PATH,
        progress_interval_mb=0,
    )
    if nbytes_rm:
        print("ReadMe downloaded")
    else:
        print("ERROR: ReadMe download failed.")

    print("\n--- STEP 3: redMaPPer DR8 cluster catalogue (VizieR TSV) ---")
    nbytes_rm_cat = download_vizier_tsv(
        "J/ApJ/785/104/table1",
        REDMAPPER_PATH,
    )
    if nbytes_rm_cat:
        print(f"redMaPPer DR8 downloaded: {nbytes_rm_cat // 1024} KB")
    else:
        print("ERROR: redMaPPer DR8 download failed.")

    print("\n--- STEP 4: redMaPPer ReadMe ---")
    nbytes_rr, _ = download_file(
        [
            "https://cdsarc.cds.unistra.fr/ftp/J/ApJ/785/104/ReadMe",
            "https://cdsarc.u-strasbg.fr/ftp/J/ApJ/785/104/ReadMe",
        ],
        REDMAPPER_README_PATH,
        progress_interval_mb=0,
    )
    if nbytes_rr:
        print("ReadMe downloaded")
    else:
        print("ERROR: redMaPPer ReadMe download failed.")

    print("\n--- STEP 5: Inspect catalogues ---")

    print("\n[Mao 2017 voids]")
    if os.path.exists(MAO_VOIDS_PATH):
        col_names = [
            "Sample", "ID", "RAdeg", "DEdeg", "z", "NGal",
            "V", "Reff", "nmin", "delmin", "r", "Prob", "Dbound",
        ]
        df_voids = pd.read_fwf(
            MAO_VOIDS_PATH,
            colspecs=[
                (0, 11), (12, 17), (18, 25), (26, 32), (33, 38),
                (39, 45), (46, 55), (56, 63), (64, 73), (74, 80),
                (81, 86), (87, 96), (97, 104),
            ],
            names=col_names,
        )
        print(f"Mao 2017 columns: {list(df_voids.columns)}")
        print(f"Total rows: {len(df_voids)}")
        print(f"First 3 rows:")
        print(df_voids.head(3).to_string())
    else:
        print("  File not found — download failed.")

    print("\n[redMaPPer DR8]")
    if os.path.exists(REDMAPPER_PATH):
        df_rm = pd.read_csv(REDMAPPER_PATH, sep="\t", comment="#", low_memory=False)
        df_rm = df_rm[~df_rm.iloc[:, 0].astype(str).str.match(r"^\s*$|^\s*-+")].reset_index(drop=True)
        df_rm = df_rm.iloc[1:].reset_index(drop=True)
        print(f"redMaPPer columns (key): {list(df_rm.columns[:15])}")
        print(f"Total clusters: {len(df_rm)}")
        print(f"First 3 rows (key columns):")
        key_cols = ["ID", "RAJ2000", "DEJ2000", "zlambda", "lambda"]
        print(df_rm[key_cols].head(3).to_string())
    else:
        print("  File not found — download failed.")

    print("\n" + "=" * 60)
    print("Data acquisition complete.")
    print("=" * 60)


if __name__ == "__main__":
    run()
