#!/usr/bin/env bash

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate gsmap_env

set -euo pipefail
cd "$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)" || exit 1

WORKDIR="/public3/ly/SMECT/true_DLPFC/results/gsmap"
GWAS_DIR="/public3/ly/SMECT/true_DLPFC/resources/gwas"
SAMPLES=(V{1..12})
DISEASES=(SCZ ADHD ANX BIP MDD)
THREADS=10

declare -A PREV=(
    [SCZ]=0.01
    [ADHD]=0.05
    [ANX]=0.10
    [BIP]=0.01
    [MDD]=0.15
)

mkdir -p ./resources
for d in gene liftover; do
    [ -d "./resources/$d" ] && continue
    src="/public3/ly/SMECT/resources/$d"
    [ -d "$src" ] || src="/home/ly/sDESE_new/resources/$d"
    cp -R "$src" "./resources/"
done

cat > fix_nan.py <<'PYEOF'
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
PYEOF

gsmap_job() {
    local sample="$1"
    local h5ad="/public3/ly/SMECT/true_DLPFC/resources/Brain_33558695_h5ad/HS_Brain_33558695_${sample}.seuratobj.Xcounts.h5ad"
    local add_latent="${WORKDIR}/${sample}/find_latent_representations/${sample}_add_latent.h5ad"

    echo "[START] gsmap ${sample}  $(date)"
    if gsmap run_find_latent_representations \
            --workdir "${WORKDIR}" --sample_name "${sample}" \
            --input_hdf5_path "${h5ad}" --annotation "spatial_region_cell_type" \
            --data_layer 'counts' \
            --feat_hidden1 256 --feat_hidden2 128 \
            --gat_hidden1 64  --gat_hidden2 64 \
            --n_neighbors 15 --weighted_adj \
            --epochs 500 --convergence_threshold 0.0001 \
            --p_drop 0.1 --gat_lr 0.001 \
            --n_comps 300 --pearson_residuals \
        && python3 fix_nan.py "${add_latent}" \
        && gsmap run_latent_to_gene \
            --workdir "${WORKDIR}" --sample_name "${sample}" \
            --annotation "spatial_region_cell_type" \
            --latent_representation 'latent_GVAE' \
            --num_neighbour 21 --num_neighbour_spatial 151; then
        echo "[DONE]  gsmap ${sample}  $(date)"
    else
        echo "[FAIL]  gsmap ${sample}  $(date)"
        echo "FAILED gsmap: ${sample}  $(date)" >> ./failed_jobs.log
    fi
}

assoc_job() {
    local d="$1" gwas="$2" sample="$3"
    local disease; disease="$(basename "${gwas}" .gwas.sum.tsv.gz)"
    local prev="${PREV[$d]:-0.01}"
    local out="./results/${d}/${disease}/${sample}"
    local score="${WORKDIR}/${sample}/latent_to_gene/${sample}_gene_marker_score.feather"
    mkdir -p "${out}"

    echo "[START] ${sample} / ${disease}  $(date)"
    if JDK_JAVA_OPTIONS=--add-opens=java.base/java.nio=ALL-UNNAMED \
        java -Xmx36g -jar /public3/ly/SMECT/kggsum.jar assoc \
            --ref-gty-file "/public3/ly/SMECT/resources/plink_gty_files/1KG.EUR.hg19.gecc.pack" refG=hg19 \
            --sum-file "${gwas}" \
              cp12Cols=chr,bp pbsCols=p betaType=1 prevalence="${prev}" refG=hg19 \
              exclude=chr6:27477797~35448354 \
            --threads 10 \
            --gene-model-database refgene \
            --output "${out}" \
            --joint-mapping method=dese \
            --gene-score-file "${score}" calcSpecificity=F \
            --gene-multiple-testing bhfdr \
            --gene-p-cut 0.05; then
        echo "[DONE]  ${sample} / ${disease}  $(date)"
    else
        echo "[FAIL]  ${sample} / ${disease}  $(date)"
        echo "FAILED: ${sample}/${disease}  $(date)" >> ./failed_jobs.log
    fi
}

fifo=".batch_fifo.$$"
mkfifo "${fifo}"
exec 6<>"${fifo}"
rm -f "${fifo}"
for ((i=0; i<THREADS; i++)); do echo; done >&6

for sample in "${SAMPLES[@]}"; do
    read -u6
    { gsmap_job "${sample}"; echo >&6; } &
done
wait

shopt -s nullglob
for d in "${DISEASES[@]}"; do
    for gwas in "${GWAS_DIR}/${d}"/*.gwas.sum.tsv.gz; do
        for sample in "${SAMPLES[@]}"; do
            read -u6
            { assoc_job "${d}" "${gwas}" "${sample}"; echo >&6; } &
        done
    done
done
wait
exec 6>&-

echo "All done.  $(date)"
