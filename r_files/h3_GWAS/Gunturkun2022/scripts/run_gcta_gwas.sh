#!/usr/bin/env bash
set -euo pipefail

ROOT="/Volumes/external_1000GB_all/playground/enhancer/r_files/h3_GWAS/Gunturkun2022"
GENO="/Volumes/external_1000GB_all/playground/enhancer/data/hs_data/v4/HS_genotypes_v4"
GCTA_DIR="/Volumes/external_1000GB_all/playground/enhancer/tools/gcta-1.95.2-macOS-arm64/gcta-1.95.2-macOS-arm64"
GCTA="$GCTA_DIR/run.sh"
PY="/Users/pete/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3"
THREADS="${THREADS:-10}"
MODE="${1:-smoke}"
TRAIT="${2:-10424}"

DATA="$ROOT/data"
RESULTS="$ROOT/results"
LOGS="$ROOT/logs"
mkdir -p "$DATA" "$RESULTS" "$LOGS"

run_gcta() {
  "$GCTA" "$@"
}

prepare_inputs() {
  if [[ ! -s "$DATA/HSNIH-PalmerPublish.csv" ]]; then
    curl -fL "https://genenetwork.org/api/v_pre1/sample_data/HSNIH-PalmerPublish.csv" \
      -o "$DATA/HSNIH-PalmerPublish.csv"
  fi

  "$PY" "$ROOT/scripts/prepare_gunturkun2022_traits.py" \
    --phenotype-csv "$DATA/HSNIH-PalmerPublish.csv" \
    --fam "$GENO.fam" \
    --out-dir "$DATA/prepared"
}

make_grm() {
  local out_prefix="$1"
  local keep_file="$2"
  local chr_arg=("${@:3}")
  run_gcta \
    --bfile "$GENO" \
    --autosome-num 20 \
    --autosome \
    ${chr_arg[@]+"${chr_arg[@]}"} \
    --maf 0.01 \
    --keep "$keep_file" \
    --make-grm \
    --out "$out_prefix" \
    --thread-num "$THREADS"
}

run_trait_mlma() {
  local trait="$1"
  local grm_prefix="$2"
  local suffix="$3"
  local mlma_mode="$4"
  local extra_args=("${@:5}")
  local pheno="$DATA/prepared/trait_${trait}.pheno"
  local out_prefix="$RESULTS/trait_${trait}${suffix}"

  run_gcta \
    --bfile "$GENO" \
    --autosome-num 20 \
    --autosome \
    --pheno "$pheno" \
    --grm "$grm_prefix" \
    "$mlma_mode" \
    ${extra_args[@]+"${extra_args[@]}"} \
    --out "$out_prefix" \
    --thread-num "$THREADS"

  local mlma_file="$out_prefix.mlma"
  if [[ ! -s "$mlma_file" && -s "$out_prefix.loco.mlma" ]]; then
    mlma_file="$out_prefix.loco.mlma"
  fi

  "$PY" "$ROOT/scripts/mlma_to_magma.py" \
    --mlma "$mlma_file" \
    --out "$out_prefix.magma.tsv"
}

case "$MODE" in
  prepare)
    prepare_inputs
    ;;
  smoke)
    prepare_inputs
    make_grm "$RESULTS/smoke_chr1_grm" "$DATA/prepared/trait_${TRAIT}.keep" --chr 1
    run_trait_mlma "$TRAIT" "$RESULTS/smoke_chr1_grm" ".smoke_chr1" --mlma --chr 1
    ;;
  full-grm)
    prepare_inputs
    make_grm "$RESULTS/target_traits_union_grm" "$DATA/prepared/target_traits_union.keep"
    ;;
  full-trait)
    prepare_inputs
    if [[ ! -s "$RESULTS/target_traits_union_grm.grm.bin" ]]; then
      make_grm "$RESULTS/target_traits_union_grm" "$DATA/prepared/target_traits_union.keep"
    fi
    run_trait_mlma "$TRAIT" "$RESULTS/target_traits_union_grm" "" --mlma-loco
    ;;
  full-all)
    prepare_inputs
    if [[ ! -s "$RESULTS/target_traits_union_grm.grm.bin" ]]; then
      make_grm "$RESULTS/target_traits_union_grm" "$DATA/prepared/target_traits_union.keep"
    fi
    for trait_id in $(seq 10418 10440); do
      if [[ -s "$RESULTS/trait_${trait_id}.magma.tsv" ]]; then
        echo "Skipping trait ${trait_id}; MAGMA input already exists."
      else
        run_trait_mlma "$trait_id" "$RESULTS/target_traits_union_grm" "" --mlma-loco
      fi
    done
    ;;
  *)
    echo "Usage: $0 {prepare|smoke [trait]|full-grm|full-trait TRAIT|full-all}" >&2
    exit 2
    ;;
esac
