#!/bin/bash

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate gsmap_env

set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR" || exit 1

mkdir -p ./resources
[ -d "./resources/gene" ] || cp -R "/public3/ly/SMECT/resources/gene" "./resources/"
[ -d "./resources/liftover" ] || cp -R "/public3/ly/SMECT/resources/liftover" "./resources/"

BASE_DIR="/public3/ly/SMECT/kggsnv_simulation"
REF_GTY_FILE="/public3/ly/SMECT/resources/plink_gty_files/ALL.hg38.maf_0.05_1k.gecc"

THREADS=10

SCENES=(
    "sim_1010_9090 25"
    "sim_1010_9090 50"
    "sim_1010_9090 75"

    "sim_3030 25"
    "sim_3030 50"
    "sim_3030 75"

    "sim_r40-59_c40-59 25"
    "sim_r40-59_c40-59 50"
    "sim_r40-59_c40-59 75"

    "sim_r20-39_c20-39_r50-79_c50-79_r90-99_c90-99 25"
    "sim_r20-39_c20-39_r50-79_c50-79_r90-99_c90-99 50"
    "sim_r20-39_c20-39_r50-79_c50-79_r90-99_c90-99 75"
)

process_file() {
    local scenario="$1"
    local ratio="$2"
    local gs_file="$3"

    local gs_basename
    gs_basename=$(basename "${gs_file}")

    local file_number=""
    if [[ "${gs_basename}" =~ ^s([0-9]+)_ ]]; then
        file_number="${BASH_REMATCH[1]}"
    else
        echo "Skip file: ${gs_file}"
        return
    fi

    local assoc_dir="${BASE_DIR}/simulateS/${scenario}/${ratio}/assoc/assoc_${file_number}"
    local sum_file="${assoc_dir}/variants.hg38.tsv.bgz"
    local output_dir="${BASE_DIR}/sdese/${scenario}/${ratio}"
    local output_prefix="${output_dir}/s${file_number}"

    mkdir -p "${output_dir}"

    JDK_JAVA_OPTIONS=--add-opens=java.base/java.nio=ALL-UNNAMED \
    java -Xmx36g -jar /public3/ly/SMECT/kggsum.jar \
            assoc \
            --ref-gty-file "${REF_GTY_FILE}" refG=hg38 \
            --sum-file "${sum_file}" \
              cp12Cols=CHROM,POS,REF,ALT \
              pbsCols=Pheno_1_Assoc@Logistic_Add_P,Pheno_1_Assoc@Logistic_Add_Beta,Pheno_1_Assoc@Logistic_Add_Beta_SE \
              betaType=1 \
              prevalence=0.01 \
              refG=hg38 \
              exclude=chr6:27477797~35448354 \
            --threads 10 \
            --gene-model-database refgene \
            --output "${output_prefix}" \
            --joint-mapping method=dese \
            --gene-score-file "${gs_file}" calcSpecificity=F \
            --gene-multiple-testing bhfdr \
            --gene-p-cut 0.05
}

for scene in "${SCENES[@]}"; do
    read -r scenario ratio <<< "${scene}"

    GSSTXT_DIR="${BASE_DIR}/sldsc/${scenario}/${ratio}/gsmap"

    tempfifo="batch_temp.$$"
    mkfifo "${tempfifo}"
    exec 6<>"${tempfifo}"
    rm -f "${tempfifo}"

    for ((i=0; i<THREADS; i++)); do
        echo
    done >&6

    for file_num in {0..99}; do
        gs_file="${GSSTXT_DIR}/s${file_num}/s${file_num}/latent_to_gene/s${file_num}_gene_marker_score.feather"

        read -u6
        {
            process_file "${scenario}" "${ratio}" "${gs_file}"
            echo >&6
        } &
    done

    wait
    exec 6>&-

done