#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="${ENHANCER_PROJECT_DIR:-/Users/pete/Desktop/playground/enhancer}"
ANALYSIS_DIR="${RESUBMIT_ANALYSIS_DIR:-${PROJECT_DIR}/r_files/revision/revision_main}"
ANALYSIS_SCRIPT="${ANALYSIS_DIR}/01_promoter_enhancer_interaction_resubmit_loop_annotation.R"
BUILDER_SCRIPT="${ANALYSIS_DIR}/build_promoter_window_anchor_sensitivity_workbook.mjs"
R_BIN="${R_BIN:-/usr/local/bin/Rscript}"
CODEX_NODE="${CODEX_NODE:-/Users/pete/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/bin/node}"
CODEX_NODE_MODULES="${CODEX_NODE_MODULES:-/Users/pete/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules}"
OUTPUT_DIR="${PROJECT_DIR}/outputs/promoter-window-anchor-sensitivity"
WORKBOOK_PATH="${OUTPUT_DIR}/promoter_window_anchor_sensitivity_atac50bp_pooled_5k_10k_25k_gene_lists.xlsx"
RUN_ROOT="$(mktemp -d /tmp/promoter_window_anchor_sensitivity.XXXXXX)"
PREVIEW_DIR="${RUN_ROOT}/previews"

cleanup() {
  rm -rf "${RUN_ROOT}"
}
trap cleanup EXIT

for required_path in \
  "${ANALYSIS_SCRIPT}" \
  "${BUILDER_SCRIPT}" \
  "${R_BIN}" \
  "${CODEX_NODE}" \
  "${CODEX_NODE_MODULES}"; do
  if [[ ! -e "${required_path}" ]]; then
    printf 'Missing required path: %s\n' "${required_path}" >&2
    exit 1
  fi
done

for width in 2001 3001 4001; do
  flank=$(( (width - 1) / 2 ))
  for mode in interval midpoint; do
    case_dir="${RUN_ROOT}/${width}_${mode}"
    mkdir -p "${case_dir}"
    printf 'Running promoter window %s bp with %s anchor matching...\n' \
      "${width}" "${mode}"
    env \
      ENHANCER_PROJECT_DIR="${PROJECT_DIR}" \
      RESUBMIT_ANALYSIS_DIR="${ANALYSIS_DIR}" \
      RESUBMIT_OUTPUT_DIR="${case_dir}" \
      PROMOTER_WINDOW_FLANK_BP="${flank}" \
      PROMOTER_ANCHOR_MATCH_MODE="${mode}" \
      ATAC_MINIMUM_OVERLAP_BP=50 \
      PROMOTER_WINDOW_SENSITIVITY_ONLY=1 \
      PROMOTER_WINDOW_GENE_LIST_ONLY=1 \
      "${R_BIN}" "${ANALYSIS_SCRIPT}"
    mv \
      "${case_dir}/promoter_window_${width}bp_anchor_${mode}_atac_overlap_50bp_gene_counts_for_workbook.tsv" \
      "${RUN_ROOT}/"
  done
done

cp "${BUILDER_SCRIPT}" "${RUN_ROOT}/build_workbook.mjs"
ln -s "${CODEX_NODE_MODULES}" "${RUN_ROOT}/node_modules"

SENSITIVITY_SOURCE_DIR="${RUN_ROOT}" \
SENSITIVITY_WORKBOOK_PATH="${WORKBOOK_PATH}" \
SENSITIVITY_PREVIEW_DIR="${PREVIEW_DIR}" \
  "${CODEX_NODE}" "${RUN_ROOT}/build_workbook.mjs"

if [[ ! -s "${WORKBOOK_PATH}" ]]; then
  printf 'Workbook was not created: %s\n' "${WORKBOOK_PATH}" >&2
  exit 1
fi

# artifact-tool may emit a large inspection sidecar next to the workbook; it is
# an internal QA artifact and not part of the requested deliverable.
rm -f "${WORKBOOK_PATH}.inspect.ndjson"

# Remove only the superseded promoter-window gene-list exports after the
# consolidated workbook has been created successfully. Comparison and QC files
# are retained.
while IFS= read -r -d '' old_file; do
  base_name="$(basename "${old_file}")"
  if [[ "${base_name}" == *_ensembl_gene_list.txt ]] ||
     [[ "${base_name}" == *_ensembl_gene_loop_counts.tsv ]] ||
     [[ "${base_name}" =~ ^promoter_window_[0-9]+bp(_anchor_(interval|midpoint))?_atac_overlap_[0-9]+bp_gene_loop_counts\.tsv$ ]]; then
    rm "${old_file}"
  fi
done < <(
  find "${ANALYSIS_DIR}" -type f \
    -path '*/sensitivity_promoter_window_*/*' \
    \( -name '*ensembl_gene_list.txt' \
       -o -name '*ensembl_gene_loop_counts.tsv' \
       -o -name 'promoter_window_*gene_loop_counts.tsv' \) \
    -print0
)

printf 'Created consolidated workbook: %s\n' "${WORKBOOK_PATH}"
