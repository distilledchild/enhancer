#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
input_dir="${script_dir}/inputs/public"
mkdir -p "${input_dir}"

download() {
  local url="$1"
  local output="$2"
  if [[ -s "${output}" ]]; then
    printf 'Using existing %s\n' "${output}"
    return
  fi
  printf 'Downloading %s\n' "${url}"
  curl -L --fail --retry 3 --retry-delay 2 --output "${output}.part" "${url}"
  mv "${output}.part" "${output}"
}

geo_root="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE193nnn/GSE193757"

download \
  "${geo_root}/soft/GSE193757_family.soft.gz" \
  "${input_dir}/GSE193757_family.soft.gz"
download \
  "${geo_root}/suppl/GSE193757_tsr.pfc.bed.gz" \
  "${input_dir}/GSE193757_tsr.pfc.bed.gz"
download \
  "${geo_root}/suppl/GSE193757_gene_counts.txt.gz" \
  "${input_dir}/GSE193757_gene_counts.txt.gz"
download \
  "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5820nnn/GSM5820551/suppl/GSM5820551_filtered_peak_bc_matrix.h5" \
  "${input_dir}/GSM5820551_filtered_peak_bc_matrix.h5"

gzip -t "${input_dir}/GSE193757_family.soft.gz"
gzip -t "${input_dir}/GSE193757_tsr.pfc.bed.gz"
gzip -t "${input_dir}/GSE193757_gene_counts.txt.gz"

(
  cd "${input_dir}"
  shasum -a 256 \
    GSE193757_family.soft.gz \
    GSE193757_tsr.pfc.bed.gz \
    GSE193757_gene_counts.txt.gz \
    GSM5820551_filtered_peak_bc_matrix.h5 \
    > input_checksums.sha256
)

printf 'Public inputs are ready in %s\n' "${input_dir}"
