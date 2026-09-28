#!/bin/bash

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate gsmap_env

set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR" || exit 1

mkdir -p ./resources
[ -d "./resources/gene" ] || cp -R "/public2/ly/SMECT/resources/gene" "./resources/"
[ -d "./resources/liftover" ] || cp -R "/public2/ly/SMECT/resources/liftover" "./resources/"

# prevalence=0.01
JDK_JAVA_OPTIONS=--add-opens=java.base/java.nio=ALL-UNNAMED \
java -Xmx96g -jar /public2/ly/SMECT/kggsum.jar \
            assoc \
            --ref-gty-file "/home/ly/sDESE_new/resources/1kg_genotype/1KG.EUR.hg19.gecc.pack" refG=hg19 \
            --sum-file "/public2/ly/SMECT/true_mouse/gwas/DBP.gwas.sum.tsv.gz" \
            cp12Cols=chr,bp,a1,a2 pbsCols=p,beta,se betaType=0 refG=hg19 exclude=chr6:27477797~35448354 \
            --threads 32 \
            --gene-model-database refgene \
            --downstream-distance 10000 \
            --upstream-distance 10000 \
            --output "/public2/ly/SMECT/true_mouse/sdese/DBP" \
            --joint-mapping method=dese \
            --gene-score-file "/public2/ly/SMECT/true_mouse/sldsc/gsmap/E16.5/latent_to_gene/E16.5_gene_marker_score.feather" \
            calcSpecificity=F \
            --gene-multiple-testing bhfdr \
            --gene-p-cut 0.05




