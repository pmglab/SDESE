import matplotlib, pandas as pd, os
matplotlib.use('Agg')

base = r"C:\Users\Administrator\Desktop\SDESE_summary\simulation"
OUT  = os.path.join(base, "metrics_precision.tsv")

SCENES = [
    "sim_1010_9090 25", "sim_1010_9090 50", "sim_1010_9090 75",
    "sim_3030 25", "sim_3030 50", "sim_3030 75",
    "sim_r20-39_c20-39_r50-79_c50-79_r90-99_c90-99 25",
    "sim_r20-39_c20-39_r50-79_c50-79_r90-99_c90-99 50",
    "sim_r20-39_c20-39_r50-79_c50-79_r90-99_c90-99 75",
    "sim_r40-59_c40-59 25", "sim_r40-59_c40-59 50", "sim_r40-59_c40-59 75",
]
SLICES = [f"s{i}" for i in range(100)]
CASES = [
    ("p-value_Condi",   "{s}_genes.hg38.condi.assoc.tsv"),
    ("spatial_Condi", "{s}_gene_marker_score.feather.genes.hg38.condi.assoc.tsv"),
]
P_COL = "Condi.ECS.P"

# "sim_1010_9090/25" -> "sim_1010_9090 25"
thr = pd.read_csv(os.path.join(base, "pvalue_thresholds.tsv"), sep="\t")
thr["scene"] = thr["scene"].str.replace("/", " ", regex=False)
thr = thr.set_index("scene").apply(pd.to_numeric, errors="coerce").fillna(0.05)

def precision(fpath, true_set, p_thr):
    try:
        df = pd.read_csv(fpath, sep="\t")
        pred = df[df[P_COL] < p_thr]["RegionID"].tolist()
        return len(set(pred) & true_set) / len(pred) if pred else None
    except Exception:
        return None

results = {}
for scene in SCENES:
    sim_dir, param = scene.split()
    true_path = f"{base}/{sim_dir}/{param}/susceptibility_genes.tsv"
    try:
        true_set = set(pd.read_csv(true_path, header=None).iloc[:, 0])
    except Exception as e:
        print(f"Skip {scene}: {e}")
        continue
    for name, tmpl in CASES:
        row = []
        for s in SLICES:
            p_thr = thr.loc[scene, s] if scene in thr.index else 0.05
            fpath = f"{base}/{sim_dir}/{param}/sdese/{tmpl.format(s=s)}"
            row.append(f"{precision(fpath, true_set, p_thr):.4f}" if precision(fpath, true_set, p_thr) is not None else "NA")
        results[f"{scene}|{name}"] = row

pd.DataFrame(results, index=SLICES).T.to_csv(OUT, sep="\t")
