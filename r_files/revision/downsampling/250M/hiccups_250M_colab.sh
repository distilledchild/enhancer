#!/bin/bash
set -euo pipefail

# Copy the seven ISAAC-generated .hic files to this Google Drive directory first.
ROOT="/content/drive/MyDrive/juicer_downsample_q30_250M"
JAR="/content/juicer_tools_1.22.01.jar"
SAMPLES=(592BB 607 D765A DA21A DA68A DBA9A DE8BA)
SEEDS=(20260731 20260737 20260732 20260736 20260734 20260735 20260733)

for i in "${!SAMPLES[@]}"; do
  SAMPLE=${SAMPLES[$i]}
  SEED=${SEEDS[$i]}
  INPUT="${ROOT}/${SAMPLE}/${SAMPLE}_inter_30_250M_seed${SEED}.hic"
  OUTPUT="${ROOT}/${SAMPLE}/${SAMPLE}_inter_30_250M_seed${SEED}.hiccups.5k10k25k"

  [[ -s "${INPUT}" ]]
  mkdir -p "${OUTPUT}"

  java -Xmx20g -jar "${JAR}" hiccups \
    -r 5000,10000,25000 \
    -k KR \
    --ignore-sparsity \
    "${INPUT}" "${OUTPUT}"

  for FILE in \
    merged_loops.bedpe \
    postprocessed_pixels_5000.bedpe \
    postprocessed_pixels_10000.bedpe \
    postprocessed_pixels_25000.bedpe; do
    [[ -s "${OUTPUT}/${FILE}" ]]
  done

  echo "${SAMPLE}: HiCCUPS 5/10/25 kb completed"
done
