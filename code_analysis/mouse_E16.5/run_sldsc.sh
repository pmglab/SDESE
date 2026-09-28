#!/bin/bash

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate gsmap_env

set -euo pipefail

WORK_DIR="/public2/ly/SMECT/true_mouse/sldsc"
GSMAP_WORKDIR="$WORK_DIR/gsmap"
SAMPLE_NAME="E16.5"
HOMOLOG_FILE="/public2/ly/SMECT/resources/mouse_human_homologs.txt"
ASSOC_FILE="/public2/ly/SMECT/true_mouse/gwas/SCZ.gwas.sum.tsv.gz"
SUMSTATS_FILE="$WORK_DIR/SCZ.sumstats.gz"
H5AD_PATH="/public2/ly/SMECT/true_mouse/E16.5_E1S1.MOSTA.h5ad"
BFILE_ROOT="/public2/ly/SMECT/resources/ref_chrom/EUR_hg19_plink/EUR.hg19.QC"
KEEP_SNP_ROOT="/public2/ly/SMECT/resources/chrom_snp/hm"
GTF_FILE="/public2/ly/SMECT/resources/gencode.v36lift37.basic.annotation.gtf"

mkdir -p "$WORK_DIR" "$GSMAP_WORKDIR"

gsmap run_find_latent_representations \
            --workdir "$GSMAP_WORKDIR" \
            --sample_name "$SAMPLE_NAME" \
            --input_hdf5_path "$H5AD_PATH" \
            --annotation annotation \
            --data_layer count

gsmap run_latent_to_gene \
            --workdir "$GSMAP_WORKDIR" \
            --sample_name "$SAMPLE_NAME" \
            --annotation annotation \
            --latent_representation latent_GVAE \
            --homolog_file "$HOMOLOG_FILE" \
            --num_neighbour 51 \
            --num_neighbour_spatial 201

zcat "$ASSOC_FILE" |
    awk '
        BEGIN {
            OFS = "\t"
            print "SNP", "A1", "A2", "Z", "N"
        }
        NR == 1 {
            for (i = 1; i <= NF; i++) {
                if ($i == "snp") snp = i
                if ($i == "a1") a1 = i
                if ($i == "a2") a2 = i
                if ($i == "z") z = i
                if ($i == "n_eff") n = i
            }
            next
        }
        NF {
            print $(snp), $(a1), $(a2), $(z), $(n)
        }
    ' |
    gzip > "$SUMSTATS_FILE"


for CHROM in {1..22}; do
    echo "[${SAMPLE_NAME}] 正在处理染色体 ${CHROM}"
    gsmap run_generate_ldscore \
        --workdir "$GSMAP_WORKDIR" \
        --sample_name "$SAMPLE_NAME" \
        --chrom "$CHROM" \
        --bfile_root "${BFILE_ROOT}" \
        --gene_window_size 10000   \
        --keep_snp_root "$KEEP_SNP_ROOT" \
        --gtf_annotation_file "$GTF_FILE"
done

gsmap run_spatial_ldsc \
        --workdir "$GSMAP_WORKDIR" \
        --sample_name "$SAMPLE_NAME" \
        --trait_name SCZ \
        --sumstats_file "$SUMSTATS_FILE" \
        --num_processes 16


