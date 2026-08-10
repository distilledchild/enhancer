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

ROOT="/lustre/isaac24/proj/UTHSC0013/panjun_work/juicer"
TARGET=140000000

SAMPLES=(74AA A2DB DA08A D765A DE8BA DA68A DBA9A DA21A 607)
EXPECTED=(141993330 145939027 213727006 284296675 302314076 378608379 380788024 382917950 419977942)
SEEDS=(20260728 20260729 20260730 20260732 20260733 20260734 20260735 20260736 20260737)

TMP=""
trap '[[ -z "${TMP}" ]] || rm -f "${TMP}"' EXIT

for i in "${!SAMPLES[@]}"; do
  SAMPLE=${SAMPLES[$i]}
  EXPECT=${EXPECTED[$i]}
  SEED=${SEEDS[$i]}

  INPUT="${ROOT}/${SAMPLE}/aligned/merged_nodups.txt"
  STATS="${ROOT}/${SAMPLE}/aligned/inter_30.txt"
  OUTDIR="${ROOT}/juicer_downsample_q30_140M/${SAMPLE}/seed${SEED}"
  OUTPUT="${OUTDIR}/${SAMPLE}_merged_nodups_q30_140M_seed${SEED}.txt"
  TMP="${OUTPUT}.part.${SLURM_JOB_ID}"

  mkdir -p "${OUTDIR}"
  rm -f "${TMP}"

  [[ -s "${INPUT}" && -s "${STATS}" ]]
  grep -q -- '-s Arima' "${STATS}"
  grep -q -- '-b ' "${STATS}"

  N=$(awk '/^Hi-C Contacts:/ {gsub(/,/, "", $3); print $3; exit}' "${STATS}")
  [[ "${N}" =~ ^[0-9]+$ && "${N}" -eq "${EXPECT}" && "${N}" -ge "${TARGET}" ]]

  if [[ -s "${OUTPUT}" ]]; then
    [[ "$(wc -l < "${OUTPUT}")" -eq "${TARGET}" ]]
    echo "${SAMPLE}: existing 140M output verified"
    continue
  fi

  echo "${SAMPLE}: downsampling ${N} to ${TARGET} contacts"

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
  [[ "${COUNT}" -eq "${TARGET}" ]]
  mv "${TMP}" "${OUTPUT}"

  printf "sample\tseed\tsource_contacts\tdownsampled_contacts\toutput\n%s\t%s\t%s\t%s\t%s\n" \
    "${SAMPLE}" "${SEED}" "${N}" "${COUNT}" "${OUTPUT}" \
    > "${OUTDIR}/downsampling_qc.tsv"

  echo "${SAMPLE}: completed"
done

echo "All remaining samples completed"
