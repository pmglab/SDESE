import re
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

def _coordinates(values, pattern):
    matches = values.astype(str).str.strip().map(re.compile(pattern).fullmatch)
    valid = matches.notna()
    coords = [(int(m.group(1)), int(m.group(2))) for m in matches[valid]]
    return np.asarray(coords, dtype=int).reshape(-1, 2), valid.to_numpy()


def _files(folder_path, pattern):
    folder = Path(folder_path)
    if not folder.is_dir():
        raise NotADirectoryError(f"Folder does not exist: {folder}")
    key = lambda p: [int(x) if x.isdigit() else x for x in re.split(r"(\d+)", p.name)]
    files = sorted(folder.glob(pattern), key=key)
    if not files:
        raise ValueError(f"No files matching {pattern!r} in {folder}")
    return files


def _load_ldsc(path, alpha):
    frame = pd.read_csv(path)
    if not {"spot", "p"}.issubset(frame.columns):
        raise ValueError(f"{path.name} must contain 'spot' and 'p' columns")
    coords, valid = _coordinates(frame["spot"], r"(\d+):(\d+)$")
    p = pd.to_numeric(frame["p"], errors="coerce").to_numpy()
    return coords, multipletests(p, alpha=alpha, method="fdr_bh")[1][valid]


def _load_dese(path, _):
    try:
        frame = pd.read_csv(path, sep="\t", dtype={"Condition": str})
    except Exception as exc:
        raise ValueError(f"Failed to read DESE result file: {path}") from exc
    if not {"Condition", "Adjusted(p)"}.issubset(frame.columns):
        raise ValueError(f"{path.name} must contain 'Condition' and 'Adjusted(p)' columns")
    if len(frame) and str(frame.iloc[-1]["Condition"]).startswith("Note:"):
        frame = frame.iloc[:-1]
    coords, valid = _coordinates(frame["Condition"], r"^C?(\d+)[:_](\d+)$")
    p = pd.to_numeric(frame.loc[valid, "Adjusted(p)"], errors="coerce").to_numpy()
    return coords, p


def _distance_to_center(coords, center):
    return np.linalg.norm(coords - np.asarray(center), axis=1)


def _distance_to_regions(coords, regions):
    if not regions:
        raise ValueError("disease_regions cannot be empty")
    x, y = coords[:, 0], coords[:, 1]
    distances = np.full(len(coords), np.inf)
    for x_min, x_max, y_min, y_max in regions:
        if x_max < x_min or y_max < y_min:
            raise ValueError("Each region requires x_min <= x_max and y_min <= y_max")
        dx = np.maximum.reduce((x_min - x, np.zeros(len(coords)), x - x_max))
        dy = np.maximum.reduce((y_min - y, np.zeros(len(coords)), y - y_max))
        distances = np.minimum(distances, np.hypot(dx, dy))
    return distances


def _membership(distances, radius, sigma):
    if radius < 0 or sigma <= 0:
        raise ValueError("radius must be non-negative and sigma must be positive")
    outside = np.maximum(distances - radius, 0)
    return np.exp(-(outside**2) / (2 * sigma**2))


def _single_sensitivity(distances, significant, radius):
    return float(np.any(significant & (distances <= radius)))


def _region_sensitivity(coords, significant, regions):
    true_spots = _distance_to_regions(coords, regions) == 0
    if not true_spots.any():
        raise ValueError("No tested spots fall inside disease_regions")
    return float(np.mean(significant[true_spots]))


def _calculate(
    files, loader, distance_fn, sensitivity_fn, radius, sigma, alpha,
):
    counts, nsr, slice_sensitivity = [], [], []
    for path in files:
        coords, p = loader(path, alpha)
        distances = distance_fn(coords)
        significant = p < alpha
        counts.append(int(significant.sum()))
        slice_sensitivity.append(sensitivity_fn(coords, distances, significant))
        if significant.any():
            m = _membership(distances, radius, sigma)[significant]
            nsr.append(float(np.sqrt(np.mean((1 - m) ** 2))))
        else:
            nsr.append(np.nan)
    return counts, float(np.mean(slice_sensitivity)), nsr, slice_sensitivity


def _summary(result):
    counts, sensitivity, nsr, _ = result
    values = np.asarray(nsr, dtype=float)
    has_nsr = np.any(~np.isnan(values))
    mean_nsr = float(np.nanmean(values)) if has_nsr else np.nan
    max_nsr = float(np.nanmax(values)) if has_nsr else np.nan
    slice_rate = float(np.mean(np.asarray(counts) > 0)) if counts else np.nan
    print(
        f"significant_spot_slice_rate={slice_rate:.6g}, sensitivity={sensitivity:.6g}, "
        f"mean_NSR_across_significant_slices={mean_nsr:.6g}, "
        f"max_NSR_across_significant_slices={max_nsr:.6g}"
    )
    return result


def _run(
    folder_path, pattern, loader, center, regions, radius, sigma, alpha,
):
    if regions is None:
        distance_fn = lambda c: _distance_to_center(c, center)
        sensitivity_fn = lambda c, d, s: _single_sensitivity(d, s, radius)
    else:
        distance_fn = lambda c: _distance_to_regions(c, regions)
        sensitivity_fn = lambda c, d, s: _region_sensitivity(c, s, regions)
    return _summary(
        _calculate(_files(folder_path, pattern), loader, distance_fn, sensitivity_fn,
                   radius, sigma, alpha)
    )


def calculate_LDSC_onediseasespot_sensitivity(
    folder_path, center, radius=0, sigma=3, sig_threshold=0.05,
):
    return _run(folder_path, "*.csv.gz", _load_ldsc, center, None,
                radius, sigma, sig_threshold)


def calculate_LDSC_region_sensitivity(
    folder_path, disease_regions, radius=0, sigma=3, sig_threshold=0.05,
):
    return _run(folder_path, "*.csv.gz", _load_ldsc, None, disease_regions,
                radius, sigma, sig_threshold)


def calculate_DESE_onediseasespot_sensitivity(
    folder_path, center, radius=0, sigma=3, sig_threshold=0.05,
):
    return _run(folder_path, "*.enrichment.tsv", _load_dese, center, None,
                radius, sigma, sig_threshold)


def calculate_DESE_region_sensitivity(
    folder_path, disease_regions, radius=0, sigma=3, sig_threshold=0.05,
):
    return _run(folder_path, "*.enrichment.tsv", _load_dese, None, disease_regions,
                radius, sigma, sig_threshold)


__all__ = [
    "calculate_LDSC_onediseasespot_sensitivity",
    "calculate_LDSC_region_sensitivity",
    "calculate_DESE_onediseasespot_sensitivity",
    "calculate_DESE_region_sensitivity",
]
