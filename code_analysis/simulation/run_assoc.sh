#!/bin/bash

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate gsmap_env

set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR" || exit 1

mkdir -p ./resources
[ -d "./resources/gene" ] || cp -R "/public3/ly/SMECT/resources/gene" "./resources/"
[ -d "./resources/liftover" ] || cp -R "/public3/ly/SMECT/resources/liftover" "./resources/"

JAVA_JAR="/public3/ly/SMECT/kggsnv.jar"
BASE_DIR="/public3/ly/SMECT"
GTY_FILE="/public3/ukb/wgs/ALL.hg38.maf_0.05.gecc.pack"
DBSNP_GECC="/public3/ly/SMECT/resources/dbsnp_157.hg38.gecc"

THREADS=10

SCENES=(
    "sim_1010_9090 25"
    "sim_1010_9090 50"
    "sim_1010_9090 75"

    "sim_3030 25"
    "sim_3030 50"
    "sim_3030 75"

    "sim_r20-39_c20-39_r50-79_c50-79_r90-99_c90-99 25"
    "sim_r20-39_c20-39_r50-79_c50-79_r90-99_c90-99 50"
    "sim_r20-39_c20-39_r50-79_c50-79_r90-99_c90-99 75"

    "sim_r40-59_c40-59 25"
    "sim_r40-59_c40-59 50"
    "sim_r40-59_c40-59 75"
)

process_file() {
    local scenario="$1"
    local ratio="$2"
    local file_number="$3"

    local phenotype_file="${BASE_DIR}/kggsnv_simulation/simulateS/${scenario}/${ratio}/ExpressionPhenotypeSimulateTask/phenotype.ped.${file_number}"
    local output_dir="${BASE_DIR}/kggsnv_simulation/simulateS/${scenario}/${ratio}/assoc/assoc_${file_number}"

    mkdir -p "${output_dir}"

    java -Xmx36g \
            -jar "${JAVA_JAR}" \
            prune \
            --input-gty-file "${GTY_FILE}" \
            --ped-file "${phenotype_file}" \
            --output "${output_dir}" \
            --threads 10 \
            --hwe 1E-5 \
            --min-obs-rate 0.8 \
            --allele-num "2~4" \
            --rsid-database dbsnp \
            "path=${DBSNP_GECC}" \
            --gene-feature-included "0~14" \
            --p-cut 1 \
            --assoc method=logistic-add
}

for scene in "${SCENES[@]}"; do
    read -r scenario ratio <<< "${scene}"

    tempfifo="batch_temp.$$"
    mkfifo "${tempfifo}"
    exec 6<>"${tempfifo}"
    rm -f "${tempfifo}"

    for ((i=0; i<THREADS; i++)); do
        echo
    done >&6

    for file_num in {0..99}; do
        phenotype_file="${BASE_DIR}/kggsnv_simulation/simulateS/${scenario}/${ratio}/ExpressionPhenotypeSimulateTask/phenotype.ped.${file_num}"

        read -u6
        {
            process_file "${scenario}" "${ratio}" "${file_num}"
            echo >&6
        } &
    done

    wait
    exec 6>&-

done