import contextlib
import io
import math
from concurrent.futures import ProcessPoolExecutor, as_completed
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd

from sensitivity_and_nsr import (
    calculate_DESE_onediseasespot_sensitivity,
    calculate_DESE_region_sensitivity,
    calculate_LDSC_onediseasespot_sensitivity,
    calculate_LDSC_region_sensitivity,
)


BASE_DIR = Path("simulation")
OUT_TSV = BASE_DIR / "sensitivity_nsr_summary.tsv"
FAILED_TSV = BASE_DIR / "sensitivity_nsr_failed.tsv"
WORKERS = 32
SIG_THRESHOLD = 0.05
FIXED_RADIUS = 3
FIXED_SIGMA = 3
DRIVER_RATIOS = [25, 50, 75]

SCENARIOS = {
    "sim_1010_9090": {
        "type": "region",
        "regions": [(10, 10, 10, 10), (90, 90, 90, 90)],
    },
    "sim_3030": {"type": "spot", "center": (30, 30)},
    "sim_r20-39_c20-39_r50-79_c50-79_r90-99_c90-99": {
        "type": "region",
        "regions": [(20, 39, 20, 39), (50, 79, 50, 79), (90, 99, 90, 99)],
    },
    "sim_r40-59_c40-59": {
        "type": "region",
        "regions": [(40, 59, 40, 59)],
    },
}

METHODS = {
    "SLDSC": {
        "subdir": "sldsc",
        "spot": calculate_LDSC_onediseasespot_sensitivity,
        "region": calculate_LDSC_region_sensitivity,
    },
    "SDESE": {
        "subdir": "sdese",
        "spot": calculate_DESE_onediseasespot_sensitivity,
        "region": calculate_DESE_region_sensitivity,
    },
}


def validate_simulation_files(folder, method):
    if method == "SLDSC":
        expected = {f"s{i}_sim.csv.gz" for i in range(100)}
        actual = {p.name for p in folder.glob("*.csv.gz")}
    else:
        expected = {
            f"s{i}_gene_marker_score.feather.enrichment.tsv" for i in range(100)
        }
        actual = {p.name for p in folder.glob("*.enrichment.tsv")}

    missing, unexpected = sorted(expected - actual), sorted(actual - expected)
    if missing or unexpected:
        raise ValueError(
            f"Simulation files do not match s0-s99; missing={missing}; "
            f"unexpected={unexpected}"
        )


def format_nsr(values):
    return ";".join("NA" if math.isnan(x) else f"{x:.10g}" for x in values)


def format_sensitivity(values):
    return ";".join(f"{x:.10g}" for x in values)


def format_counts(values):
    return ";".join(str(int(x)) for x in values)


def summarize_nsr(values):
    values = np.asarray(values, dtype=float)
    if np.all(np.isnan(values)):
        return np.nan, np.nan
    return float(np.nanmean(values)), float(np.nanmax(values))


def failed_row(method, scenario, ratio, folder, error):
    return {
        "failed": True,
        "method": method,
        "scenario": scenario,
        "driver_ratio": ratio,
        "folder": str(folder),
        "error": error,
    }


def run_one(task):
    scenario, ratio, method = task
    scenario_info, method_info = SCENARIOS[scenario], METHODS[method]
    folder = BASE_DIR / scenario / str(ratio) / method_info["subdir"]
    if not folder.is_dir():
        return failed_row(method, scenario, ratio, folder, "Result folder does not exist")

    try:
        validate_simulation_files(folder, method)
        kwargs = {
            "folder_path": folder,
            "radius": FIXED_RADIUS,
            "sigma": FIXED_SIGMA,
            "sig_threshold": SIG_THRESHOLD,
        }
        key = "center" if scenario_info["type"] == "spot" else "disease_regions"
        kwargs[key] = scenario_info[key if key == "center" else "regions"]
        func = method_info[scenario_info["type"]]

        with contextlib.redirect_stdout(io.StringIO()):
            counts, sensitivity, nsr, simulation_sensitivity = func(**kwargs)

        counts = np.asarray(counts, dtype=int)
        nsr = np.asarray(nsr, dtype=float)
        simulation_sensitivity = np.asarray(simulation_sensitivity, dtype=float)
        for name, values in {
            "significant-spot counts": counts,
            "NSR values": nsr,
            "simulation sensitivity values": simulation_sensitivity,
        }.items():
            if len(values) != 100:
                raise ValueError(f"Expected 100 {name}, but received {len(values)}")
    except Exception as exc:
        error = f"{type(exc).__name__}: {exc}"
        return failed_row(method, scenario, ratio, folder, error)

    mean_nsr, max_nsr = summarize_nsr(nsr)
    return {
        "method": method,
        "scenario": scenario,
        "driver_ratio": ratio,
        "mean_sensitivity": sensitivity,
        "mean_NSR": mean_nsr,
        "max_NSR": max_nsr,
        "simulation_significant_spot_counts": format_counts(counts),
        "simulation_NSR_values": format_nsr(nsr),
        "simulation_sensitivity_values": format_sensitivity(simulation_sensitivity),
        "significant_spot_simulation_rate": float(np.mean(counts > 0)),
    }


def main():
    tasks = list(product(SCENARIOS, DRIVER_RATIOS, METHODS))
    rows, failed = [], []
    with ProcessPoolExecutor(max_workers=WORKERS) as pool:
        futures = [pool.submit(run_one, task) for task in tasks]
        for i, future in enumerate(as_completed(futures), 1):
            try:
                row = future.result()
            except Exception as exc:
                failed.append({"error": f"{type(exc).__name__}: {exc}"})
                continue
            (failed if row.get("failed") else rows).append(row)
            if i % 10 == 0 or i == len(tasks):
                print(f"Finished {i}/{len(tasks)}")

    columns = [
        "method", "scenario", "driver_ratio", "mean_sensitivity", "mean_NSR",
        "max_NSR", "simulation_significant_spot_counts", "simulation_NSR_values",
        "simulation_sensitivity_values", "significant_spot_simulation_rate",
    ]
    result = pd.DataFrame(rows, columns=columns)
    if not result.empty:
        result = result.sort_values(["scenario", "driver_ratio", "method"])
    result.to_csv(OUT_TSV, sep="\t", index=False)
    print(f"Saved {len(result)} rows to {OUT_TSV}")

    if failed:
        pd.DataFrame(failed).to_csv(FAILED_TSV, sep="\t", index=False)
        print(f"Failed {len(failed)} tasks; details saved to {FAILED_TSV}")


if __name__ == "__main__":
    main()
