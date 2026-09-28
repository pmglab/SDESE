import sys, anndata as ad, numpy as np, scipy.sparse as sp
p = sys.argv[1]
a = ad.read_h5ad(p)
pr = a.layers['pearson_residuals']
arr = pr.toarray() if sp.issparse(pr) else np.array(pr)
n_nan = int(np.isnan(arr).sum())
n_inf = int(np.isinf(arr).sum())
if n_nan or n_inf:
    arr = np.nan_to_num(arr, nan=0.0, posinf=0.0, neginf=0.0)
    a.layers['pearson_residuals'] = arr
    a.write_h5ad(p)
    print(f"[fix] {p}: {n_nan} NaN, {n_inf} Inf -> 0")
