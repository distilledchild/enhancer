#!/bin/bash
#SBATCH -A ACF-UTHSC0013
#SBATCH -p long-bigmem
#SBATCH -q long-bigmem
#SBATCH --time 23:59:59
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=64G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pkim11@uthsc.edu
#SBATCH -o /lustre/isaac24/proj/UTHSC0013/panjun_work/juicer/create_hic_140M_%j.log
#SBATCH -J hic-ds-140M

set -euo pipefail
export LC_ALL=C

ROOT="/lustre/isaac24/proj/UTHSC0013/panjun_work/juicer"
CHROM_SIZES="/lustre/isaac24/proj/UTHSC0013/panjun_work/refs/rn7_ucsc/bwa_rn7chr/chrom.sizes"
JUICER_JAR="/lustre/isaac24/proj/UTHSC0013/panjun_work/tools/juicer/CPU/common/juicer_tools.jar"
RESTRICTION_SITES="/lustre/isaac24/proj/UTHSC0013/panjun_work/refs/rn7_ucsc/juicer_restriction_sites/rn7chr_juicer_arima4.txt"
TARGET=140000000

SAMPLES=(592BB 607 74AA A2DB D765A DA08A DA21A DA68A DBA9A DE8BA)
SEEDS=(20260731 20260737 20260728 20260729 20260732 20260730 20260736 20260734 20260735 20260733)

TMP=""
trap '[[ -z "${TMP}" ]] || rm -f "${TMP}"' EXIT

for FILE in "${JUICER_JAR}" "${RESTRICTION_SITES}" "${CHROM_SIZES}"; do
  [[ -s "${FILE}" ]] || { echo "Missing shared input: ${FILE}" >&2; exit 1; }
done

for i in "${!SAMPLES[@]}"; do
  SAMPLE=${SAMPLES[$i]}
  SEED=${SEEDS[$i]}
  OUTDIR="${ROOT}/juicer_downsample_q30_140M/${SAMPLE}/seed${SEED}"
  INPUT="${OUTDIR}/${SAMPLE}_merged_nodups_q30_140M_seed${SEED}.txt"
  OUTPUT="${OUTDIR}/${SAMPLE}_inter_30_140M_seed${SEED}.hic"
  TMP="${OUTPUT}.part.${SLURM_JOB_ID}"

  [[ -s "${INPUT}" ]] || { echo "Missing input: ${INPUT}" >&2; exit 1; }
  [[ "$(wc -l < "${INPUT}")" -eq "${TARGET}" ]]

  if [[ -s "${OUTPUT}" ]]; then
    [[ "$(head -c 3 "${OUTPUT}")" == "HIC" && "$(stat -c '%s' "${OUTPUT}")" -gt 1000000 ]]
    echo "${SAMPLE}: existing .hic verified"
    continue
  fi

  rm -f "${TMP}"
  java -Djava.awt.headless=true -Xms8g -Xmx56g \
    -jar "${JUICER_JAR}" pre -q 0 -f "${RESTRICTION_SITES}" \
    -r 5000,10000,25000 "${INPUT}" "${TMP}" "${CHROM_SIZES}"

  [[ -s "${TMP}" && "$(head -c 3 "${TMP}")" == "HIC" && "$(stat -c '%s' "${TMP}")" -gt 1000000 ]]
  mv "${TMP}" "${OUTPUT}"
  printf "sample\tcontacts\tresolutions\thic_file\n%s\t%s\t%s\t%s\n" \
    "${SAMPLE}" "${TARGET}" "5kb,10kb,25kb" "${OUTPUT}" \
    > "${OUTDIR}/hic_creation_qc.tsv"
done
