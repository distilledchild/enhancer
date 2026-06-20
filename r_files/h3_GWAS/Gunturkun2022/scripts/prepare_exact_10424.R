#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: prepare_exact_10424.R COHORT_TSV FAM OUT_DIR", call. = FALSE)
}

cohort_tsv <- args[[1]]
fam_path <- args[[2]]
out_dir <- args[[3]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

rint <- function(x) {
  qnorm((rank(x, ties.method = "average") - 0.5) / length(x))
}

clean_num <- function(x) {
  x[x %in% c("", "x", "X", "NA", "NaN", "nan", "-9")] <- NA
  suppressWarnings(as.numeric(x))
}

cohort <- read.delim(cohort_tsv, check.names = FALSE, stringsAsFactors = FALSE)
names(cohort) <- sub("^([0-9]+):.*$", "HSR_\\1", names(cohort))
names(cohort)[names(cohort) == "sample_id"] <- "id"

fam <- read.table(fam_path, stringsAsFactors = FALSE)
fam_ids <- fam[[2]]
cohort <- cohort[cohort$id %in% fam_ids, ]

y_raw <- clean_num(cohort$HSR_10424)
sex <- clean_num(cohort$HSR_10443)
batch <- clean_num(cohort$HSR_10444)
color <- clean_num(cohort$HSR_10445)
center <- clean_num(cohort$HSR_10446)

ok <- complete.cases(y_raw, sex, batch, color, center)
dat <- data.frame(
  id = cohort$id[ok],
  y_raw = y_raw[ok],
  sex = factor(sex[ok]),
  batch = factor(batch[ok]),
  color = factor(color[ok]),
  center = factor(center[ok]),
  stringsAsFactors = FALSE
)

dat$y_sex_rint <- NA_real_
for (level in levels(dat$sex)) {
  idx <- dat$sex == level
  dat$y_sex_rint[idx] <- rint(dat$y_raw[idx])
}

covariates <- c("batch", "color")
covar_rows <- list()
selected <- c()
for (covar in covariates) {
  fit <- lm(as.formula(paste("y_sex_rint ~", covar)), data = dat)
  fit0 <- lm(y_sex_rint ~ 1, data = dat)
  an <- anova(fit0, fit)
  p_value <- an$`Pr(>F)`[2]
  r2 <- summary(fit)$r.squared
  include <- is.finite(p_value) && p_value < 0.05 && r2 > 0.02
  if (include) {
    selected <- c(selected, covar)
  }
  covar_rows[[covar]] <- data.frame(
    covariate = covar,
    levels = length(levels(dat[[covar]])),
    p_value = p_value,
    r2 = r2,
    selected = include
  )
}

if (length(selected) > 0) {
  formula <- as.formula(paste("y_sex_rint ~", paste(selected, collapse = " + ")))
  fit <- lm(formula, data = dat)
  dat$y_resid <- resid(fit) + mean(dat$y_sex_rint)
} else {
  dat$y_resid <- dat$y_sex_rint
}
dat$y_final <- rint(dat$y_resid)

pheno <- data.frame(FID = dat$id, IID = dat$id, PHENO = signif(dat$y_final, 12))
keep <- data.frame(FID = dat$id, IID = dat$id)

write.table(pheno, file.path(out_dir, "trait_10424.exact_public.pheno"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(keep, file.path(out_dir, "trait_10424.exact_public.keep"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(dat, file.path(out_dir, "trait_10424.exact_public.preprocess_values.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(do.call(rbind, covar_rows),
            file.path(out_dir, "trait_10424.exact_public.covariate_tests.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

summary_rows <- data.frame(
  metric = c(
    "trait",
    "n_after_fam_and_complete_covariates",
    "female_0",
    "male_1",
    "selected_covariates",
    "missing_age_covariate",
    "phenotype_file",
    "keep_file"
  ),
  value = c(
    "10424 open_field_totaldistance",
    nrow(dat),
    sum(dat$sex == "0"),
    sum(dat$sex == "1"),
    ifelse(length(selected) > 0, paste(selected, collapse = "+"), "none"),
    "true",
    file.path(out_dir, "trait_10424.exact_public.pheno"),
    file.path(out_dir, "trait_10424.exact_public.keep")
  )
)
write.table(summary_rows, file.path(out_dir, "trait_10424.exact_public.preprocess_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("wrote exact-public 10424 phenotype\n")
cat("n:", nrow(dat), "\n")
cat("selected covariates:", ifelse(length(selected) > 0, paste(selected, collapse = "+"), "none"), "\n")
