import os
import gzip
import numpy as np
import pandas as pd
import requests
from astropy.io import fits

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")


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


def run():
    print("=" * 60)
    print("VOID-DOMAIN v1 — Real Data Acquisition")
    print("=" * 60)

    print("\n--- STEP 1: Mao et al. (2017) BOSS DR12 void catalogue ---")
    voids_path = os.path.join(DATA_DIR, "mao2017_boss_voids.dat")
    voids_urls = [
        "https://cdsarc.cds.unistra.fr/ftp/J/ApJ/835/161/table1.dat",
        "https://cdsarc.u-strasbg.fr/ftp/J/ApJ/835/161/table1.dat",
    ]
    nbytes, used_url = download_file(voids_urls, voids_path, progress_interval_mb=0)
    if nbytes:
        print(f"Mao 2017 void catalogue downloaded: {nbytes // 1024} KB")
    else:
        print("ERROR: All Mao 2017 void catalogue URLs failed.")

    print("\n--- STEP 2: VizieR ReadMe ---")
    readme_path = os.path.join(DATA_DIR, "mao2017_ReadMe.txt")
    readme_urls = [
        "https://cdsarc.cds.unistra.fr/ftp/J/ApJ/835/161/ReadMe",
        "https://cdsarc.u-strasbg.fr/ftp/J/ApJ/835/161/ReadMe",
    ]
    nbytes_rm, _ = download_file(readme_urls, readme_path, progress_interval_mb=0)
    if nbytes_rm:
        print("ReadMe downloaded")
    else:
        print("ERROR: ReadMe download failed.")

    print("\n--- STEP 3: redMaPPer DR8 cluster catalogue ---")
    rm_path_v63 = os.path.join(DATA_DIR, "redmapper_dr8_v6.3.fits.gz")
    rm_path_v510 = os.path.join(DATA_DIR, "redmapper_dr8_v5.10.fits.gz")

    nbytes_rm_cat, used_rm_url = download_file(
        [
            "http://risa.stanford.edu/redmapper/redmapper_dr8_public_v6.3_catalog.fits.gz",
            "https://risa.stanford.edu/redmapper/redmapper_dr8_public_v6.3_catalog.fits.gz",
        ],
        rm_path_v63,
        progress_interval_mb=10,
    )
    rm_path = rm_path_v63
    if not nbytes_rm_cat:
        print("  v6.3 failed, trying v5.10...")
        nbytes_rm_cat, used_rm_url = download_file(
            [
                "http://risa.stanford.edu/redmapper/redmapper_dr8_public_v5.10_catalog.fits.gz",
                "https://risa.stanford.edu/redmapper/redmapper_dr8_public_v5.10_catalog.fits.gz",
            ],
            rm_path_v510,
            progress_interval_mb=10,
        )
        rm_path = rm_path_v510

    if nbytes_rm_cat:
        print(f"redMaPPer DR8 downloaded: {nbytes_rm_cat / (1024 * 1024):.1f} MB")
    else:
        print("ERROR: All redMaPPer DR8 URLs failed.")
        rm_path = None

    print("\n--- STEP 4: Inspect catalogues ---")

    print("\n[Mao 2017 voids]")
    if os.path.exists(voids_path):
        try:
            df_voids = pd.read_csv(voids_path, sep=r"\s+")
            print(f"Mao 2017 columns: {list(df_voids.columns)}")
            print(f"Total rows: {len(df_voids)}")
            print(f"First 3 rows:")
            print(df_voids.head(3).to_string())
        except Exception as e:
            print(f"  Parse error: {e}")
            print("  Trying fixed-width format...")
            try:
                with open(voids_path) as f:
                    first_lines = [f.readline() for _ in range(5)]
                print("  Raw first 5 lines:")
                for line in first_lines:
                    print(f"    {line.rstrip()}")
            except Exception as e2:
                print(f"  Raw read also failed: {e2}")
    else:
        print("  File not found — download failed.")

    print("\n[redMaPPer DR8]")
    if rm_path and os.path.exists(rm_path):
        try:
            hdulist = fits.open(rm_path)
            cols = hdulist[1].columns.names
            n_clusters = len(hdulist[1].data)
            print(f"redMaPPer columns: {cols}")
            print(f"Total clusters: {n_clusters}")
            df_rm = pd.DataFrame({c: hdulist[1].data[c] for c in cols[:10]})
            print(f"First 3 rows (first 10 columns):")
            print(df_rm.head(3).to_string())
            hdulist.close()
        except Exception as e:
            print(f"  Parse error: {e}")
    else:
        print("  File not found — download failed.")

    print("\n" + "=" * 60)
    print("Data acquisition complete.")
    print("=" * 60)


if __name__ == "__main__":
    run()
