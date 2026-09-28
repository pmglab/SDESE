#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR" || exit 1

declare -A HOT_LOCATIONS=(
  [sim_1010_9090]="10:10,90:90"
  [sim_3030]="30:30"
  [sim_r20-39_c20-39_r50-79_c50-79_r90-99_c90-99]="20~39:20~39,50~79:50~79,90~99:90~99"
  [sim_r40-59_c40-59]="40~59:40~59"
)

declare -A MAX_DRIVER_NUM=(
  [25]="125"
  [50]="250"
  [75]="375"
)

SCENARIOS=(
  sim_1010_9090
  sim_3030
  sim_r20-39_c20-39_r50-79_c50-79_r90-99_c90-99
  sim_r40-59_c40-59
)

PROPORTIONS=(25 50 75)

run_simulation() {
  local scenario="$1"
  local proportion="$2"
  local output_dir="${SCRIPT_DIR}/simulateS/${scenario}/${proportion}"

  mkdir -p "${output_dir}"
  echo "Running ${scenario}/${proportion}: output=${output_dir}"

  java -Xmx56g -jar "/public3/ly/SMECT/kggsnv.jar" \
    simulate \
      --input-gty-file "/public3/ukb/wgs/ALL.hg38.maf_0.05.gecc.pack" refG="hg38" \
      --xqtl-gene-file "/public3/ly/SMECT/resources/eqtlgenes.txt" \
      --output "${output_dir}" \
      --threads "24" \
      --min-obs-rate "0.4" \
      --allele-num "2~4" \
      --local-maf "0.05~0.5" \
      --r2-cut "0.01" \
      --gene-feature-included "0~14" \
      --xqtl-gene maxQTLProportion="0.4" \
                  effectScale="1" \
                  heritability="0.2" \
      --phenotype-gene maxDriverNum="${MAX_DRIVER_NUM[${proportion}]}" \
                       heritability="0.15" \
                       prevalence="0.01" \
      --sampling-group groupNum="100" \
      --sampling sampleSize="8000" caseProportion="0.5" \
      --spatial-gene chipSize="100" \
                     hotLocations="${HOT_LOCATIONS[${scenario}]}" \
                     maxHotExpressionGeneNum="500" \
                     dropoutRate="0.4" \
                     chipNum="100"
}

for scenario in "${SCENARIOS[@]}"; do
  for proportion in "${PROPORTIONS[@]}"; do
    run_simulation "${scenario}" "${proportion}"
  done
done
