import sys

print("VOID-DOMAIN v1 — ready")

if __name__ == "__main__":
    if "--mock-calibration" in sys.argv:
        from void_domain.mock_calibration import (
            download_tng300_data,
            build_mock_void_table,
            run_mock_fit,
            finalise_thresholds,
        )

        print("\n========== MOCK CALIBRATION MODE ==========\n")

        download_tng300_data()
        df = build_mock_void_table()
        results = run_mock_fit(df)
        finalise_thresholds(results)
