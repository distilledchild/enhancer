#!/usr/bin/env bash
set -euo pipefail

PROJECT="/Volumes/external_1000GB_all/playground/enhancer/data/Khun_2025/gwas_inputs/u01_peter_kalivas"
GENO="/Volumes/external_1000GB_all/playground/enhancer/data/Khun_2025/hs_rats_round10_2"
TOOLS="/Volumes/external_1000GB_all/playground/enhancer/tools/kuhn_repro"
PLINK="$TOOLS/plink-1.90b6.21/bin/plink"
GCTA_DIR="$TOOLS/gcta-1.94.2-macOS-arm64/gcta-1.94.2-MacOS-ARM-x86_64"
GCTA="$GCTA_DIR/gcta64"
export DYLD_LIBRARY_PATH="$GCTA_DIR/gcclib:${DYLD_LIBRARY_PATH:-}"

THREADS="${THREADS:-4}"
TARGET_FILTERED_SNPS="${TARGET_FILTERED_SNPS:-5424892}"
KEEP_IDS="$PROJECT/genotypes/keep_rfids.txt"
KEEP_PLINK="$PROJECT/genotypes/keep_fid_iid_plink.txt"
SEX_UPDATE="$PROJECT/genotypes/update_sex_plink.txt"

mkdir -p "$PROJECT"/{logs,genotypes,grm,results/gwas,temp}

timestamp() {
  date "+%Y-%m-%d %H:%M:%S"
}

log() {
  printf "[%s] %s\n" "$(timestamp)" "$*"
}

require_file() {
  if [[ ! -e "$1" ]]; then
    echo "Missing required file: $1" >&2
    exit 1
  fi
}

print_versions() {
  log "Tool versions"
  "$PLINK" --version
  { "$GCTA" 2>&1 || true; } | sed -n '1,7p'
}

run_plink_qc() {
  local name="$1"
  shift
  local out="$PROJECT/genotypes/$name"

  if [[ -s "$out.lmiss" && -s "$out.hwe" && -s "$out.frq" && "${FORCE:-0}" != "1" ]]; then
    log "Skipping existing PLINK QC: $name"
    return
  fi

  log "Running PLINK QC: $name"
  "$PLINK" \
    --bfile "$GENO" \
    "$@" \
    --hardy \
    --keep "$KEEP_PLINK" \
    --update-sex "$SEX_UPDATE" \
    --thread-num "$THREADS" \
    --freq \
    --missing \
    --nonfounders \
    --out "$out" \
    --chr-set 20 no-xy
}

make_accepted_snps() {
  local accepted="$PROJECT/genotypes/accepted_snps.txt"
  local quality="$PROJECT/genotypes/snpquality.tsv"
  local thresholds="$PROJECT/genotypes/parameter_thresholds.txt"

  if [[ -s "$accepted" && "${FORCE:-0}" != "1" ]]; then
    local existing_n
    existing_n="$(wc -l < "$accepted" | tr -d ' ')"
    if [[ "$existing_n" == "$TARGET_FILTERED_SNPS" ]]; then
      log "Skipping existing accepted SNP list: $existing_n matches report target"
      return
    fi
    log "Existing accepted SNP list has $existing_n SNPs; rebuilding"
  fi

  log "Combining PLINK QC outputs and writing accepted SNP list"
  : > "$accepted"
  printf "SNP\tCHR\tF_MISS\tHWE\tMAF\tPASS_MISS\tPASS_MAF\tPASS_HWE\tPASS\n" > "$quality"

  for name in autosomes xfilter yfilter; do
    paste \
      <(tail -n +2 "$PROJECT/genotypes/$name.lmiss") \
      <(tail -n +2 "$PROJECT/genotypes/$name.hwe") \
      <(tail -n +2 "$PROJECT/genotypes/$name.frq") |
    awk -v miss_thr=0.1 -v hwe_thr=1e-10 -v maf_thr=0.005 -v accepted="$accepted" -v quality="$quality" '
      BEGIN { OFS="\t" }
      {
        chr=$1
        snp=$2
        fmiss=$5 + 0
        hwe=$14 + 0
        maf=$19 + 0

        pass_miss = ((fmiss < miss_thr) || (chr == 22))
        pass_maf = ((maf >= maf_thr) && (maf <= (1 - maf_thr)))
        pass_hwe = ((hwe > hwe_thr) || (chr == 22) || (chr == 24))
        pass = (pass_miss && pass_maf && pass_hwe)

        print snp, chr, fmiss, $14, maf, pass_miss, pass_maf, pass_hwe, pass >> quality
        if (pass) {
          print snp >> accepted
        }
      }
    '
  done

  gzip -f "$quality"
  {
    echo "--geno 0.1"
    echo "--maf 0.005"
    echo "--hwe 1e-10"
  } > "$thresholds"

  local n
  n="$(wc -l < "$accepted" | tr -d ' ')"
  log "Accepted SNP count: $n; report target: $TARGET_FILTERED_SNPS"
}

make_filtered_bed() {
  local out="$PROJECT/genotypes/genotypes"

  if [[ -s "$out.bed" && -s "$out.bim" && -s "$out.fam" && "${FORCE:-0}" != "1" ]]; then
    log "Skipping existing filtered PLINK dataset"
    return
  fi

  log "Creating filtered PLINK dataset"
  "$PLINK" \
    --bfile "$GENO" \
    --extract "$PROJECT/genotypes/accepted_snps.txt" \
    --keep "$KEEP_PLINK" \
    --update-sex "$SEX_UPDATE" \
    --make-bed \
    --thread-num "$THREADS" \
    --set-missing-var-ids '@:#' \
    --keep-allele-order \
    --set-hh-missing \
    --make-founders \
    --out "$out" \
    --chr-set 20 no-xy

  local n
  n="$(wc -l < "$out.bim" | tr -d ' ')"
  log "Filtered BIM SNP count: $n; report target: $TARGET_FILTERED_SNPS"
}

filter_stage() {
  require_file "$PLINK"
  require_file "$GCTA"
  require_file "$GENO.bed"
  require_file "$GENO.bim"
  require_file "$GENO.fam"
  require_file "$KEEP_IDS"

  print_versions
  log "Using THREADS=$THREADS"
  prepare_plink_sample_files
  log "Keep sample count: $(wc -l < "$KEEP_PLINK" | tr -d ' ')"

  run_plink_qc autosomes --chr 1-20 MT
  run_plink_qc xfilter --chr X --filter-females
  run_plink_qc yfilter --chr Y --filter-males
  make_accepted_snps
  make_filtered_bed
}

prepare_plink_sample_files() {
  log "Preparing PLINK FID/IID keep file and sex update file"
  python3 - <<'PY'
import csv
from pathlib import Path

project = Path("/Volumes/external_1000GB_all/playground/enhancer/data/Khun_2025/gwas_inputs/u01_peter_kalivas")
geno = Path("/Volumes/external_1000GB_all/playground/enhancer/data/Khun_2025/hs_rats_round10_2")
processed = Path("/Volumes/external_1000GB_all/playground/enhancer/data/Khun_2025/gwas_inputs/processed_data_ready.csv")

keep_ids = project / "genotypes" / "keep_rfids.txt"
keep_plink = project / "genotypes" / "keep_fid_iid_plink.txt"
sex_update = project / "genotypes" / "update_sex_plink.txt"

iid_to_fid = {}
with open(str(geno) + ".fam") as f:
    for line in f:
        fields = line.split()
        if len(fields) >= 2:
            iid_to_fid[fields[1]] = fields[0]

sex_by_iid = {}
with open(processed, newline="") as f:
    reader = csv.DictReader(f)
    for row in reader:
        iid = row.get("rfid", "")
        sex = row.get("sex", "")
        if sex in {"M", "m", "male", "1"}:
            sex_by_iid[iid] = "1"
        elif sex in {"F", "f", "female", "2"}:
            sex_by_iid[iid] = "2"

kept = []
with open(keep_ids) as f:
    for line in f:
        parts = line.split()
        if not parts:
            continue
        iid = parts[-1]
        fid = iid_to_fid.get(iid)
        if fid is None:
            raise SystemExit(f"RFID not found in FAM IID column: {iid}")
        kept.append((fid, iid))

with open(keep_plink, "w") as out_keep, open(sex_update, "w") as out_sex:
    for fid, iid in kept:
        out_keep.write(f"{fid} {iid}\n")
        sex = sex_by_iid.get(iid)
        if sex is not None:
            out_sex.write(f"{fid} {iid} {sex}\n")

print(f"wrote {len(kept)} keep rows to {keep_plink}")
print(f"wrote {sum(1 for _ in open(sex_update))} sex rows to {sex_update}")
PY
}

grm_stage() {
  local geno="$PROJECT/genotypes/genotypes"
  require_file "$GCTA"
  require_file "$geno.bed"
  require_file "$geno.bim"
  require_file "$geno.fam"

  print_versions
  log "Building per-chromosome GRMs with THREADS=$THREADS"
  : > "$PROJECT/grm/listofchrgrms.txt"

  for chr in $(seq 1 20); do
    local out="$PROJECT/grm/${chr}chrGRM"
    if [[ -s "$out.grm.bin" && -s "$out.grm.N.bin" && -s "$out.grm.id" && "${FORCE:-0}" != "1" ]]; then
      log "Skipping existing chr$chr GRM"
    else
      log "Making chr$chr GRM"
      "$GCTA" \
        --thread-num "$THREADS" \
        --bfile "$geno" \
        --chr "$chr" \
        --autosome-num 20 \
        --make-grm-bin \
        --out "$out"
    fi
    echo "$out" >> "$PROJECT/grm/listofchrgrms.txt"
  done

  log "Merging autosomal GRMs"
  "$GCTA" \
    --thread-num "$THREADS" \
    --mgrm "$PROJECT/grm/listofchrgrms.txt" \
    --make-grm-bin \
    --out "$PROJECT/grm/AllchrGRM"
}

gwas_one_stage() {
  local trait="${1:-regressedlr_total_heroin_consumption}"
  local geno="$PROJECT/genotypes/genotypes"
  local pheno="$PROJECT/data/pheno/$trait.txt"

  require_file "$GCTA"
  require_file "$geno.bed"
  require_file "$PROJECT/grm/AllchrGRM.grm.bin"
  require_file "$pheno"

  print_versions
  log "Running Kuhn-style per-chromosome MLMA for trait: $trait"
  for chr in $(seq 1 20); do
    local out="$PROJECT/results/gwas/${trait}_chrgwas${chr}"
    if [[ -s "$out.mlma" && "${FORCE:-0}" != "1" ]]; then
      log "Skipping existing GWAS chr$chr for $trait"
      continue
    fi
    log "GWAS chr$chr for $trait"
    "$GCTA" \
      --thread-num 1 \
      --pheno "$pheno" \
      --bfile "$geno" \
      --grm "$PROJECT/grm/AllchrGRM" \
      --autosome-num 20 \
      --chr "$chr" \
      --mlma-subtract-grm "$PROJECT/grm/${chr}chrGRM" \
      --mlma \
      --out "$out"
  done
}

case "${1:-filter}" in
  filter)
    filter_stage
    ;;
  grm)
    grm_stage
    ;;
  gwas-one)
    gwas_one_stage "${2:-regressedlr_total_heroin_consumption}"
    ;;
  all-one)
    filter_stage
    grm_stage
    gwas_one_stage "${2:-regressedlr_total_heroin_consumption}"
    ;;
  *)
    echo "Usage: $0 {filter|grm|gwas-one [trait]|all-one [trait]}" >&2
    exit 2
    ;;
esac
