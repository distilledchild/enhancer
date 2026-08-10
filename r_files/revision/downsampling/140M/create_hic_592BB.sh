#!/bin/bash
#SBATCH -A ACF-UTHSC0013
#SBATCH -p long-bigmem
#SBATCH -q long-bigmem
#SBATCH --time 23:59:59
#SBATCH -N 1
#SBATCH -n 30
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pkim11@uthsc.edu
#SBATCH -o /lustre/isaac24/proj/UTHSC0013/panjun_work/juicer/create_hic_%x_%j.log
#SBATCH -J hic-ds-592BB-1.6

eval "$(conda shell.bash hook)"

set -euo pipefail
export LC_ALL=C

SAMPLE="592BB"
SEED=20260731
TARGET=140000000

ROOT="/lustre/isaac24/proj/UTHSC0013/panjun_work/juicer"
OUTDIR="${ROOT}/juicer_downsample_q30_140M/${SAMPLE}/seed${SEED}"
INPUT="${OUTDIR}/${SAMPLE}_merged_nodups_q30_140M_seed${SEED}.txt"
OUTPUT="${OUTDIR}/${SAMPLE}_inter_30_140M_seed${SEED}.hic"
TMP="${OUTPUT}.part.${SLURM_JOB_ID}"

JUICER_JAR="${ROOT}/${SAMPLE}/juicer_tools.jar"
RESTRICTION_SITES="${ROOT}/${SAMPLE}/restriction_sites/rn7chr_juicer_arima4.txt"
CHROM_SIZES="/lustre/isaac24/proj/UTHSC0013/panjun_work/refs/rn7_ucsc/bwa_rn7chr/chrom.sizes"

trap 'rm -f "${TMP}"' EXIT

for file in "${INPUT}" "${JUICER_JAR}" "${RESTRICTION_SITES}" "${CHROM_SIZES}"; do
  [[ -s "${file}" ]] || {
    echo "Missing input: ${file}" >&2
    exit 1
  }
done

COUNT=$(wc -l < "${INPUT}")
[[ "${COUNT}" -eq "${TARGET}" ]] || {
  echo "Expected ${TARGET} contacts but found ${COUNT}" >&2
  exit 1
}

if [[ -s "${OUTPUT}" ]]; then
  echo "Output already exists: ${OUTPUT}"
  exit 0
fi

java -Djava.awt.headless=true -Xms8g -Xmx56g \
  -jar "${JUICER_JAR}" pre \
  -q 0 \
  -f "${RESTRICTION_SITES}" \
  -r 5000,10000,25000 \
  "${INPUT}" "${TMP}" "${CHROM_SIZES}"

[[ -s "${TMP}" ]] || {
  echo "Juicer Tools did not create a .hic file" >&2
  exit 1
}

[[ "$(head -c 3 "${TMP}")" == "HIC" && "$(stat -c '%s' "${TMP}")" -gt 1000000 ]] || {
  echo "The generated file failed the .hic format/size check" >&2
  exit 1
}

mv "${TMP}" "${OUTPUT}"

printf "sample\tcontacts\tresolutions\thic_file\n%s\t%s\t%s\t%s\n" \
  "${SAMPLE}" "${TARGET}" "5kb,10kb,25kb" "${OUTPUT}" \
  > "${OUTDIR}/hic_creation_qc.tsv"

echo "Completed: ${OUTPUT}"

