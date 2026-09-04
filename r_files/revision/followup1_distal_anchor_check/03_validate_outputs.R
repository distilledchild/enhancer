# Validate internal consistency and fixed phase-1 analysis invariants.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

current_script_path <- function() {
  file.args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (!length(file.args)) return(normalizePath(getwd(), mustWork = TRUE))
  normalizePath(sub("^--file=", "", file.args[[1]]), mustWork = TRUE)
}

script.dir <- dirname(current_script_path())
result.dir <- file.path(script.dir, "results")

required.files <- c(
  "distal_anchor_evidence.tsv.gz",
  "loop_gene_evidence.tsv.gz",
  "complete_four_layer_loop_gene_candidates.tsv.gz",
  "gene_evidence.tsv.gz",
  "evidence_summary_by_resolution.tsv",
  "matched_null_resamples.tsv.gz",
  "matched_null_summary.tsv",
  "matched_control_availability.tsv",
  "expression_mapping_summary.tsv",
  "expression_mapping_exclusions.tsv",
  "expression_association_test.tsv",
  "hic_support_balance.tsv",
  "pooled_pfc_tsr_summary.tsv",
  "snatac_rn6_to_rn7_audit.tsv",
  "analysis_parameters.tsv"
)
missing.files <- required.files[!file.exists(file.path(result.dir, required.files))]
if (length(missing.files)) {
  stop("Missing required result file(s): ", paste(missing.files, collapse = ", "))
}

distal <- read_tsv(
  file.path(result.dir, "distal_anchor_evidence.tsv.gz"),
  show_col_types = FALSE
)
loop.gene <- read_tsv(
  file.path(result.dir, "loop_gene_evidence.tsv.gz"),
  show_col_types = FALSE
)
complete.candidates <- read_tsv(
  file.path(result.dir, "complete_four_layer_loop_gene_candidates.tsv.gz"),
  show_col_types = FALSE
)
summary.by.resolution <- read_tsv(
  file.path(result.dir, "evidence_summary_by_resolution.tsv"),
  show_col_types = FALSE
)
null.resamples <- read_tsv(
  file.path(result.dir, "matched_null_resamples.tsv.gz"),
  show_col_types = FALSE
)
null.summary <- read_tsv(
  file.path(result.dir, "matched_null_summary.tsv"),
  show_col_types = FALSE
)
match.availability <- read_tsv(
  file.path(result.dir, "matched_control_availability.tsv"),
  show_col_types = FALSE
)
expression.mapping <- read_tsv(
  file.path(result.dir, "expression_mapping_summary.tsv"),
  show_col_types = FALSE
)
tsr.summary <- read_tsv(
  file.path(result.dir, "pooled_pfc_tsr_summary.tsv"),
  show_col_types = FALSE
)
snatac.audit <- read_tsv(
  file.path(result.dir, "snatac_rn6_to_rn7_audit.tsv"),
  show_col_types = FALSE
)
parameters <- read_tsv(
  file.path(result.dir, "analysis_parameters.tsv"),
  show_col_types = FALSE
)

lookup_value <- function(df, key.column, value.column, key) {
  values <- df[[value.column]][df[[key.column]] == key]
  if (length(values) != 1L) return(NA)
  values[[1]]
}

checks <- list()
add_check <- function(name, observed, expected, pass) {
  checks[[length(checks) + 1L]] <<- tibble(
    check = name,
    observed = as.character(observed),
    expected = as.character(expected),
    pass = isTRUE(pass)
  )
}

core.columns <- c(
  "loop_id", "resolution", "anchor_chr", "anchor_start", "anchor_end",
  "atac_overlap_ge50", "pooled_pfc_distal_tsr_gt3kb_present",
  "colocated_atac_distal_tsr", "complete_four_layer_support"
)

add_check(
  "primary loop rows", nrow(distal), 12295L,
  nrow(distal) == 12295L
)
add_check(
  "unique primary loop IDs", n_distinct(distal$loop_id), nrow(distal),
  n_distinct(distal$loop_id) == nrow(distal)
)
add_check(
  "core columns present", sum(core.columns %in% names(distal)), length(core.columns),
  all(core.columns %in% names(distal))
)
if (all(core.columns %in% names(distal))) {
  add_check(
    "missing core values", sum(is.na(distal[, core.columns])), 0L,
    !anyNA(distal[, core.columns])
  )
}
add_check(
  "ATAC-supported loop count", sum(distal$atac_overlap_ge50), 10469L,
  sum(distal$atac_overlap_ge50) == 10469L
)
add_check(
  "colocated evidence requires same-anchor ATAC and distal TSR",
  sum(distal$colocated_atac_distal_tsr & !distal$atac_and_distal_tsr_same_anchor),
  0L,
  !any(distal$colocated_atac_distal_tsr & !distal$atac_and_distal_tsr_same_anchor)
)
add_check(
  "unique loop-gene pairs", n_distinct(paste(loop.gene$loop_id, loop.gene$gene_id)),
  nrow(loop.gene),
  n_distinct(paste(loop.gene$loop_id, loop.gene$gene_id)) == nrow(loop.gene)
)
add_check(
  "complete-candidate loop count",
  n_distinct(complete.candidates$loop_id), sum(distal$complete_four_layer_support),
  n_distinct(complete.candidates$loop_id) == sum(distal$complete_four_layer_support)
)
add_check(
  "complete-candidate evidence flags",
  all(
    complete.candidates$colocated_atac_distal_tsr &
      complete.candidates$naive_pfc_expression_support
  ),
  TRUE,
  all(
    complete.candidates$colocated_atac_distal_tsr &
      complete.candidates$naive_pfc_expression_support
  )
)

expected.resolutions <- c("5K", "10K", "25K", "ALL")
add_check(
  "resolution summary rows", nrow(summary.by.resolution), 4L,
  nrow(summary.by.resolution) == 4L
)
add_check(
  "resolution labels", paste(summary.by.resolution$resolution, collapse = ","),
  paste(expected.resolutions, collapse = ","),
  identical(summary.by.resolution$resolution, expected.resolutions)
)
all.summary <- summary.by.resolution %>% filter(resolution == "ALL")
add_check(
  "ALL summary loop count", all.summary$n_one_sided_promoter_anchored_loops,
  nrow(distal),
  nrow(all.summary) == 1L &&
    all.summary$n_one_sided_promoter_anchored_loops == nrow(distal)
)

expected.metrics <- c(
  "atac_overlap_ge50",
  "pooled_pfc_tsr_present",
  "pooled_pfc_distal_tsr_gt3kb_present",
  "atac_and_distal_tsr_same_anchor",
  "colocated_atac_distal_tsr"
)
expected.null.rows <- length(expected.resolutions) * length(expected.metrics)
add_check(
  "matched-null summary combinations", nrow(null.summary), expected.null.rows,
  nrow(null.summary) == expected.null.rows &&
    n_distinct(paste(null.summary$resolution, null.summary$metric)) == expected.null.rows
)
configured.resamples <- as.integer(lookup_value(
  parameters, "parameter", "value", "n_resamples"
))
resamples.per.group <- null.resamples %>%
  count(resolution, metric, name = "n")
add_check(
  "resamples per null combination",
  paste(sort(unique(resamples.per.group$n)), collapse = ","),
  configured.resamples,
  nrow(resamples.per.group) == expected.null.rows &&
    all(resamples.per.group$n == configured.resamples)
)
all.match.availability <- match.availability %>% filter(resolution == "ALL")
add_check(
  "control-match availability accounting",
  sum(all.match.availability$n_observed_anchors), nrow(distal),
  sum(all.match.availability$n_observed_anchors) == nrow(distal) &&
    !any(all.match.availability$match_level == "no_control_available")
)
rate.columns <- c(
  "observed_rate", "mean_null_rate", "null_rate_q025", "null_rate_q975",
  "empirical_p", "fdr_bh"
)
rates.in.range <- all(vapply(
  null.summary[, rate.columns],
  function(x) all(is.finite(x) & x >= 0 & x <= 1),
  logical(1)
))
add_check("null rates and P values in [0,1]", rates.in.range, TRUE, rates.in.range)

tsr.source.total <- as.integer(lookup_value(
  tsr.summary, "metric", "n", "pooled_pfc_tsr_source_total"
))
tsr.analyzed <- as.integer(lookup_value(
  tsr.summary, "metric", "n", "pooled_pfc_tsr_analyzed"
))
tsr.excluded <- as.integer(lookup_value(
  tsr.summary, "metric", "n", "pooled_pfc_tsr_excluded_nonanalysis_seqlevels"
))
add_check(
  "TSR source accounting", tsr.analyzed + tsr.excluded, tsr.source.total,
  tsr.analyzed + tsr.excluded == tsr.source.total
)

source.features <- as.integer(lookup_value(
  snatac.audit, "metric", "value", "source_rn6_unique_features"
))
retained.features <- as.integer(lookup_value(
  snatac.audit, "metric", "value", "source_features_retained"
))
not.retained.features <- as.integer(lookup_value(
  snatac.audit, "metric", "value", "source_features_not_retained"
))
unexpected.local <- as.integer(lookup_value(
  snatac.audit, "metric", "value", "local_names_not_in_source_h5"
))
add_check(
  "snATAC source-feature accounting",
  retained.features + not.retained.features, source.features,
  retained.features + not.retained.features == source.features
)
add_check(
  "unexpected lifted snATAC source names", unexpected.local, 0L,
  unexpected.local == 0L
)
add_check(
  "exact current-symbol expression mappings do not exceed genes",
  expression.mapping$n_exact_current_gene_id_symbol_expression_matches,
  paste0("<=", expression.mapping$n_unique_ensembl_gene_ids),
  expression.mapping$n_exact_current_gene_id_symbol_expression_matches <=
    expression.mapping$n_unique_ensembl_gene_ids
)

validation <- bind_rows(checks)
write_tsv(validation, file.path(result.dir, "validation_checks.tsv"))
print(validation, n = Inf)

if (any(!validation$pass)) {
  stop("Output validation failed; inspect results/validation_checks.tsv.")
}
message("All ", nrow(validation), " output validation checks passed.")
