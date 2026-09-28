from __future__ import annotations

from pathlib import Path

import h5py
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "pdf.fonttype": 42,
    "font.size": 7,
})

ROOT = Path(__file__).resolve().parent
H5AD_PATH = ROOT / "E16.5_E1S1.MOSTA.h5ad"
OUTPUT_DIR = ROOT / "figures"


def decode(values: np.ndarray) -> list[str]:
    return [
        value.decode("utf-8") if isinstance(value, (bytes, np.bytes_)) else str(value)
        for value in values
    ]


def load_annotation() -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str], list[str]]:
    with h5py.File(H5AD_PATH, "r") as handle:
        spatial = np.asarray(handle["obsm/spatial"][:], dtype=float)
        codes = np.asarray(handle["obs/annotation/codes"][:], dtype=int)
        categories = decode(handle["obs/annotation/categories"][:])
        colors = decode(handle["uns/annotation_colors"][:])

    if spatial.shape != (len(codes), 2):
        raise ValueError(f"Unexpected spatial array shape: {spatial.shape}")
    if len(categories) != len(colors):
        raise ValueError("Tissue categories and colors have different lengths")
    if codes.min() < 0 or codes.max() >= len(categories):
        raise ValueError("Invalid tissue annotation codes")
    return spatial[:, 0], spatial[:, 1], codes, categories, colors


def make_figure() -> plt.Figure:
    x, y, codes, categories, colors = load_annotation()
    fig = plt.figure(figsize=(5.2, 5.0), dpi=150, facecolor="white")
    grid = fig.add_gridspec(1, 2, width_ratios=[1.00, 0.68], wspace=0.03)
    ax = fig.add_subplot(grid[0, 0])
    legend_ax = fig.add_subplot(grid[0, 1])

    # Draw category-by-category to keep the legend order identical to the H5AD.
    for code, color in enumerate(colors):
        mask = codes == code
        ax.scatter(
            x[mask], y[mask], s=1.15, c=color, alpha=0.92,
            linewidths=0, edgecolors="none", rasterized=True,
        )

    ax.set_aspect("equal", adjustable="box")
    ax.axis("off")
    ax.margins(0.01)

    handles = [
        Line2D([0], [0], marker="o", linestyle="none", markerfacecolor=color,
               markeredgecolor="none", markersize=3.4, label=category)
        for category, color in zip(categories, colors)
    ]
    legend_ax.axis("off")
    legend_ax.legend(
        handles=handles, loc="center left", frameon=False, fontsize=8.2,
        handletextpad=0.65, labelspacing=0.42, borderaxespad=0,
    )
    fig.subplots_adjust(left=0.01, right=0.995, top=0.995, bottom=0.01)
    return fig


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig = make_figure()
    stem = OUTPUT_DIR / "E16.5_tissue_annotation_spatial_map"
    fig.savefig(Path(f"{stem}.png"), dpi=600, bbox_inches="tight", pad_inches=0.03, facecolor="white")
    fig.savefig(Path(f"{stem}.pdf"), bbox_inches="tight", pad_inches=0.03, facecolor="white")
    plt.close(fig)
    print(f"Saved tissue annotation map for 121,767 spots and 25 tissues to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
