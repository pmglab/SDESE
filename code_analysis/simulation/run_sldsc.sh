#!/usr/bin/env bash

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate gsmap_env_gpu

set -euo pipefail

SMECT_DIR="/public3/ly/SMECT"
SIMULATION_DIR="${SMECT_DIR}/kggsnv_simulation/simulateS"
BFILE_ROOT="${SMECT_DIR}/resources/plink_gty_files/with_rsid/split_chrom/ALL.hg38.maf_0.01_1k_map_rsid.chr"
KEEP_SNP_ROOT="${SMECT_DIR}/resources/chrom_snp/hm"
GTF_FILE="${SMECT_DIR}/resources/gencode.v46.basic.annotation.gtf"

OUTPUT_ROOT_DIR="${SMECT_DIR}/kggsnv_simulation/sldsc"

SAMPLE_START=0
SAMPLE_END=99

MAX_JOBS=10
GWAS_N=8000
SLDSC_NUM_PROCESSES=12

SCENES=(
    "sim_1010_9090|25"
    "sim_1010_9090|50"
    "sim_1010_9090|75"

    "sim_3030|25"
    "sim_3030|50"
    "sim_3030|75"

    "sim_r20-39_c20-39_r50-79_c50-79_r90-99_c90-99|25"
    "sim_r20-39_c20-39_r50-79_c50-79_r90-99_c90-99|50"
    "sim_r20-39_c20-39_r50-79_c50-79_r90-99_c90-99|75"

    "sim_r40-59_c40-59|25"
    "sim_r40-59_c40-59|50"
    "sim_r40-59_c40-59|75"
)

run_parallel_samples() {
    local job_func="$1"
    local failed=0
    local pids=()

    for SAMPLE_NUM in $(seq "$SAMPLE_START" "$SAMPLE_END"); do
        (
            "$job_func" "$SAMPLE_NUM"
        ) &

        pids+=("$!")

        if ((${#pids[@]} == MAX_JOBS)); then
            failed=0
            for pid in "${pids[@]}"; do
                wait "$pid" || failed=1
            done
            ((failed == 0)) || return 1
            pids=()
        fi
    done

    failed=0
    for pid in "${pids[@]+"${pids[@]}"}"; do
        wait "$pid" || failed=1
    done

    ((failed == 0)) || return 1
}

run_scene() {
    local SCENARIO="$1"
    local SPE="$2"

    local OUTPUT_DIR="${OUTPUT_ROOT_DIR}/${SCENARIO}/${SPE}"
    local H5AD_DIR="${SIMULATION_DIR}/${SCENARIO}/${SPE}/h5ad"
    local ASSOC_DIR="${SIMULATION_DIR}/${SCENARIO}/${SPE}/assoc"
    local GSMAP_WORKDIR="${OUTPUT_DIR}/gsmap"
    local SUMSTATS_DIR="${OUTPUT_DIR}/sumstats"

    mkdir -p "$GSMAP_WORKDIR" "$SUMSTATS_DIR"

    echo "开始场景: ${SCENARIO}/${SPE}"
    echo "H5AD: ${H5AD_DIR}"
    echo "OUTPUT: ${OUTPUT_DIR}"

    step1_find_latent() {
        local SAMPLE_NUM="$1"
        local SAMPLE_NAME="s${SAMPLE_NUM}"
        local H5AD_FILE="${H5AD_DIR}/gene.spatial.expression.txt.${SAMPLE_NUM}.h5ad"

        gsmap run_find_latent_representations \
            --workdir "$GSMAP_WORKDIR/s${SAMPLE_NUM}" \
            --sample_name "$SAMPLE_NAME" \
            --input_hdf5_path "$H5AD_FILE" \
            --annotation annotation \
            --data_layer counts

    }

    step2_latent_to_gene() {
        local SAMPLE_NUM="$1"
        local SAMPLE_NAME="s${SAMPLE_NUM}"

        gsmap run_latent_to_gene \
            --workdir "$GSMAP_WORKDIR/s${SAMPLE_NUM}" \
            --sample_name "$SAMPLE_NAME" \
            --annotation annotation \
            --latent_representation latent_GVAE \
            --num_neighbour 51 \
            --num_neighbour_spatial 201

    }

    step3_make_sumstats() {
        local SAMPLE_NUM="$1"
        local SAMPLE_NAME="s${SAMPLE_NUM}"
        local ASSOC_FILE="${ASSOC_DIR}/assoc_${SAMPLE_NUM}/variants.hg38.tsv.bgz"
        local SUMSTATS_FILE="${SUMSTATS_DIR}/variants.ldsc.anno.sumstats_${SAMPLE_NUM}.gz"

        zcat "$ASSOC_FILE" |
            awk -v n="$GWAS_N" '
                BEGIN {
                    OFS = "\t"
                    print "SNP", "A1", "A2", "Z", "N"
                }
                $1 == "#CHROM" {
                    for (i = 1; i <= NF; i++) {
                        if ($i == "REF") ref = i
                        if ($i == "ALT") alt = i
                        if ($i == "rsid") snp = i
                        if ($i == "Pheno_1_Assoc@Logistic_Add_Beta") beta = i
                        if ($i == "Pheno_1_Assoc@Logistic_Add_Beta_SE") se = i
                    }
                    if (!beta || !se) {
                        print "错误: 未找到 beta/SE 列" > "/dev/stderr"
                        exit 1
                    }
                    next
                }
                !/^#/ && NF {
                    z = ($(se) == 0) ? "NA" : $(beta) / $(se)
                    print $(snp), $(alt), $(ref), z, n
                }
            ' |
            gzip > "$SUMSTATS_FILE"

    }

    step4_generate_ldscore() {
        local SAMPLE_NUM="$1"
        local SAMPLE_NAME="s${SAMPLE_NUM}"

        for CHROM in {1..22}; do
            echo "[${SCENARIO}/${SPE}/${SAMPLE_NAME}] 处理染色体 ${CHROM}"

            gsmap run_generate_ldscore \
                --workdir "$GSMAP_WORKDIR/s${SAMPLE_NUM}" \
                --sample_name "$SAMPLE_NAME" \
                --chrom "$CHROM" \
                --bfile_root "$BFILE_ROOT" \
                --keep_snp_root "$KEEP_SNP_ROOT" \
                --gtf_annotation_file "$GTF_FILE"
        done

    }

    step5_spatial_ldsc() {
        local SAMPLE_NUM="$1"
        local SAMPLE_NAME="s${SAMPLE_NUM}"
        local SUMSTATS_FILE="${SUMSTATS_DIR}/variants.ldsc.anno.sumstats_${SAMPLE_NUM}.gz"

        gsmap run_spatial_ldsc \
            --workdir "$GSMAP_WORKDIR/s${SAMPLE_NUM}" \
            --sample_name "$SAMPLE_NAME" \
            --trait_name sim \
            --sumstats_file "$SUMSTATS_FILE" \
            --num_processes "$SLDSC_NUM_PROCESSES"

    }

    run_parallel_samples step1_find_latent

    run_parallel_samples step2_latent_to_gene

    run_parallel_samples step3_make_sumstats

    run_parallel_samples step4_generate_ldscore

    run_parallel_samples step5_spatial_ldsc
}

for scene in "${SCENES[@]}"; do
    IFS='|' read -r SCENARIO SPE <<< "$scene"
    run_scene "$SCENARIO" "$SPE"
done

echo "Time: $((ELAPSED_SECONDS / 60)) m $((ELAPSED_SECONDS % 60)) s"
