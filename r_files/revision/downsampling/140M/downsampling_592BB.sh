#!/bin/bash
#SBATCH -A ACF-UTHSC0013
#SBATCH -p campus-bigmem
#SBATCH -q campus-bigmem
#SBATCH --time 23:59:59
#SBATCH -N 1
#SBATCH -n 15
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pkim11@uthsc.edu
#SBATCH -o /lustre/isaac24/proj/UTHSC0013/panjun_work/juicer/downsampling_remaining_%j.log
#SBATCH -J DS_remaining


set -euo pipefail
export LC_ALL=C

SAMPLE="592BB"
EXPECTED=259559222
TARGET=140000000
SEED=20260731

ROOT="/lustre/isaac24/proj/UTHSC0013/panjun_work/juicer"
INPUT="${ROOT}/${SAMPLE}/aligned/merged_nodups.txt"
STATS="${ROOT}/${SAMPLE}/aligned/inter_30.txt"
OUTDIR="${ROOT}/juicer_downsample_q30_140M/${SAMPLE}/seed${SEED}"
OUTPUT="${OUTDIR}/${SAMPLE}_merged_nodups_q30_140M_seed${SEED}.txt"
TMP="${OUTPUT}.part.${SLURM_JOB_ID}"

mkdir -p "${OUTDIR}"
trap 'rm -f "${TMP}"' EXIT

[[ -s "${INPUT}" && -s "${STATS}" ]] || {
  echo "Missing input for ${SAMPLE}" >&2
  exit 1
}

# Refuse to use a source run that does not record the Arima -s/-b settings.
grep -q -- '-s Arima' "${STATS}"
grep -q -- '-b ' "${STATS}"

N=$(
  awk '/^Hi-C Contacts:/ {gsub(/,/, "", $3); print $3; exit}' "${STATS}"
)

[[ "${N}" =~ ^[0-9]+$ && "${N}" -eq "${EXPECTED}" && "${N}" -ge "${TARGET}" ]] || {
  echo "Unexpected contact count for ${SAMPLE}: ${N:-missing}" >&2
  exit 1
}

if [[ -s "${OUTPUT}" ]]; then
  [[ "$(wc -l < "${OUTPUT}")" -eq "${TARGET}" ]] || {
    echo "Existing output has an incorrect line count: ${OUTPUT}" >&2
    exit 1
  }
  echo "Verified output already exists: ${OUTPUT}"
  exit 0
fi

# Exact one-pass sampling of MAPQ >=30, non-intra-fragment contacts.
gawk -v N="${N}" -v K="${TARGET}" -v seed="${SEED}" '
BEGIN {srand(seed); remaining=N; need=K; eligible=0; selected=0; bad=0; extra=0}
NF < 12 {bad++; next}
$9 >= 30 && $12 >= 30 && !($2 == $6 && $4 == $8) {
  eligible++
  if (remaining <= 0) {extra++; next}
  if (need == remaining || rand() < need / remaining) {print; need--; selected++}
  remaining--
}
END {
  if (bad || extra || eligible != N || selected != K || need != 0 || remaining != 0) exit 10
}' "${INPUT}" > "${TMP}"

COUNT=$(wc -l < "${TMP}")
[[ "${COUNT}" -eq "${TARGET}" ]] || {
  echo "Generated ${COUNT} contacts instead of ${TARGET}" >&2
  exit 1
}

mv "${TMP}" "${OUTPUT}"
printf "sample\tseed\tsource_contacts\tdownsampled_contacts\toutput\n%s\t%s\t%s\t%s\t%s\n" \
  "${SAMPLE}" "${SEED}" "${N}" "${COUNT}" "${OUTPUT}" \
  > "${OUTDIR}/downsampling_qc.tsv"

echo "Completed: ${OUTPUT}"


