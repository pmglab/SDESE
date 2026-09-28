from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd
import scanpy as sc
from statsmodels.stats.multitest import multipletests


DATA_ROOT = Path(r"C:\Users\Administrator\Desktop\SDESE_summary\mouse_E16.5")
SDESE_TRAIT = "TP"
SLDSC_TRAIT = "MDD"
ALPHA = 0.05

COLOR_MAP = LinearSegmentedColormap.from_list(
    "fire_red",
    ["#FFD700", "#FF8C00", "#FF4500", "#B22222", "#8B0000"],
    N=256,
)


def add_spatial_coordinates(adata: sc.AnnData, results: pd.DataFrame) -> pd.DataFrame:
    results = results.loc[results.index.intersection(adata.obs_names)].copy()
    obs = adata.obs.join(results, how="left")
    obs["x"] = adata.obsm["spatial"][:, 0]
    obs["y"] = adata.obsm["spatial"][:, 1]
    return obs


def split_by_significance(
    obs: pd.DataFrame,
    adjusted_p_column: str,
    alpha: float = ALPHA,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    mask = obs[adjusted_p_column].lt(alpha).fillna(False)
    significant = obs.loc[mask].copy()
    significant["neg_log_p"] = -np.log10(
        significant[adjusted_p_column].clip(lower=np.finfo(float).tiny)
    )
    significant = significant.sort_values("neg_log_p")
    non_significant = obs.loc[~mask].copy()
    return significant, non_significant


def plot_spatial_results(
    significant: pd.DataFrame,
    non_significant: pd.DataFrame,
    vmax: float,
    output_path: Path,
) -> plt.Figure:
    background_color = "#F0F8FF"
    fig, ax = plt.subplots(figsize=(3, 5), dpi=150, facecolor=background_color)
    ax.axis("off")

    ax.scatter(
        non_significant["x"],
        non_significant["y"],
        color="#6495ED",
        alpha=0.03,
        s=5,
        edgecolors="none",
        rasterized=True,
    )

    scatter = ax.scatter(
        significant["x"],
        significant["y"],
        c=significant["neg_log_p"],
        cmap=COLOR_MAP,
        vmin=1.0,
        vmax=vmax,
        s=5,
        alpha=0.85,
        edgecolors="#333333",
        linewidth=0.3,
        rasterized=True,
    )

    colorbar_axis = fig.add_axes([0.88, 0.15, 0.03, 0.7])
    colorbar = fig.colorbar(scatter, cax=colorbar_axis)
    colorbar.set_label(r"$-\log_{10}(\mathrm{adjusted}\ p)$", fontsize=11, labelpad=10)
    colorbar.ax.tick_params(labelsize=9)
    colorbar.set_ticks(np.linspace(1.0, vmax, 6 if vmax == 2.0 else 5))

    ax.text(
        0,
        1,
        f"{len(significant)} significant features",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=6,
        bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "none", "pad": 3},
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        output_path,
        format="tiff",
        dpi=600,
        bbox_inches="tight",
        pil_kwargs={"compression": "tiff_lzw"},
    )
    return fig


def run_sdese(data_root: Path, trait: str) -> plt.Figure:
    adata = sc.read_h5ad(data_root / "E16.5_E1S1.MOSTA.h5ad")
    result_path = data_root / "sdese" / (
        f"{trait}_E16.5_gene_marker_score.feather.enrichment.tsv"
    )
    results = pd.read_csv(result_path, sep="\t").set_index("Condition")
    obs = add_spatial_coordinates(adata, results)
    significant, non_significant = split_by_significance(obs, "Adjusted(p)")

    return plot_spatial_results(
        significant,
        non_significant,
        vmax=2.0,
        output_path=data_root / "figures" / f"{trait}_SDESE_spatial_visualization.tiff",
    )


def run_sldsc(data_root: Path, trait: str) -> plt.Figure:
    adata = sc.read_h5ad(data_root / "E16.5_E1S1.MOSTA.h5ad")
    results = pd.read_csv(data_root / "sldsc" / f"E16.5_{trait}.csv.gz")
    results = results.set_index("spot")

    valid_p = results["p"].notna()
    results["adjusted_p"] = np.nan
    results.loc[valid_p, "adjusted_p"] = multipletests(
        results.loc[valid_p, "p"], method="fdr_bh"
    )[1]

    obs = add_spatial_coordinates(adata, results)
    significant, non_significant = split_by_significance(obs, "adjusted_p")

    return plot_spatial_results(
        significant,
        non_significant,
        vmax=5.0,
        output_path=data_root / "figures" / f"{trait}_SLDSC_spatial_visualization_bh.tiff",
    )


# Run both analyses
run_sdese(DATA_ROOT, SDESE_TRAIT)
run_sldsc(DATA_ROOT, SLDSC_TRAIT)
plt.show()
