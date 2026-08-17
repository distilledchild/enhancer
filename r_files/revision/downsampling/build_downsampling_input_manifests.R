# Build symmetric inventories for the organized 140M and 250M Google Drive inputs.

current_script_path <- function() {
  file.args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file.args) == 0L) return(NA_character_)
  file.path <- gsub("~\\+~", " ", sub("^--file=", "", file.args[[1]]))
  normalizePath(
    file.path, winslash = "/", mustWork = FALSE
  )
}

script.path <- current_script_path()
project.root <- if (is.na(script.path)) {
  NA_character_
} else {
  dirname(dirname(dirname(dirname(script.path))))
}
bundled.root <- if (!is.na(project.root) && dir.exists(project.root)) {
  file.path(project.root, "data", "juicer_downsample_q30_140M_250M")
} else {
  NA_character_
}
configured.root <- Sys.getenv("DOWNSAMPLING_COMBINED_ROOT", unset = "")
legacy.root <- paste0(
  "~/Library/CloudStorage/GoogleDrive-wellclouder@gmail.com/My Drive/",
  "juicer_downsample_q30_140M_250M"
)
root.candidates <- path.expand(c(bundled.root, configured.root, legacy.root))
root.candidates <- root.candidates[
  !is.na(root.candidates) & nzchar(root.candidates) & dir.exists(root.candidates)
]
if (length(root.candidates) == 0L) {
  stop("Combined downsampling root is unavailable.")
}
root <- normalizePath(root.candidates[[1]], winslash = "/", mustWork = TRUE)

metadata <- data.frame(
  sample = c(
    "592BB", "607", "74AA", "A2DB", "D765A",
    "DA08A", "DA21A", "DA68A", "DBA9A", "DE8BA"
  ),
  source_contacts = c(
    259559222, 419977942, 141993330, 145939027, 284296675,
    213727006, 382917950, 378608379, 380788024, 302314076
  ),
  seed = c(
    20260731, 20260737, 20260728, 20260729, 20260732,
    20260730, 20260736, 20260734, 20260735, 20260733
  )
)

nonempty <- function(path) file.exists(path) && file.info(path)$size > 0
size_or_na <- function(path) if (file.exists(path)) file.info(path)$size else NA_real_
md5_or_na <- function(path) {
  if (!nonempty(path)) return(NA_character_)
  unname(tools::md5sum(path))
}
allocated_blocks_or_na <- function(path) {
  if (!file.exists(path)) return(NA_real_)
  stat.args <- if (identical(Sys.info()[["sysname"]], "Darwin")) {
    c("-f", "%b", shQuote(path))
  } else {
    c("-c", "%b", shQuote(path))
  }
  value <- suppressWarnings(system2(
    "stat", stat.args, stdout = TRUE, stderr = FALSE
  ))
  if (length(value) == 0L) NA_real_ else suppressWarnings(as.numeric(value[[1]]))
}

build_manifest <- function(depth) {
  target <- as.numeric(sub("M$", "", depth)) * 1000000
  rows <- lapply(seq_len(nrow(metadata)), function(i) {
    sample <- metadata$sample[[i]]
    seed <- metadata$seed[[i]]
    eligible <- metadata$source_contacts[[i]] >= target
    sample.dir <- file.path(root, depth, sample)
    prefix <- sprintf("%s_inter_30_%s_seed%s", sample, depth, seed)
    hic <- file.path(sample.dir, paste0(prefix, ".hic"))
    downsampling.qc <- file.path(sample.dir, "downsampling_qc.tsv")
    hic.qc <- file.path(sample.dir, "hic_creation_qc.tsv")
    hiccups.dir <- file.path(sample.dir, paste0(prefix, ".hiccups.5k10k25k"))
    merged <- file.path(hiccups.dir, "merged_loops.bedpe")
    hic.blocks <- if (eligible) allocated_blocks_or_na(hic) else NA_real_
    postprocessed <- file.path(
      hiccups.dir,
      paste0("postprocessed_pixels_", c(5000, 10000, 25000), ".bedpe")
    )
    complete <- eligible && nonempty(hic) && nonempty(downsampling.qc) &&
      nonempty(hic.qc) && nonempty(merged) && all(vapply(postprocessed, nonempty, logical(1)))
    data.frame(
      depth = depth,
      sample = sample,
      source_contacts = metadata$source_contacts[[i]],
      target_contacts = target,
      sampling_fraction = if (eligible) target / metadata$source_contacts[[i]] else NA_real_,
      seed = seed,
      eligible = eligible,
      status = if (complete) {
        "complete_hic_and_hiccups_5k10k25k"
      } else if (!eligible) {
        "excluded_source_contacts_below_target"
      } else {
        "eligible_but_incomplete"
      },
      hic_file = if (eligible) hic else NA_character_,
      hic_size_bytes = if (eligible) size_or_na(hic) else NA_real_,
      hic_allocated_blocks = hic.blocks,
      hic_storage_state = if (!eligible) {
        NA_character_
      } else if (is.na(hic.blocks)) {
        "unknown"
      } else if (hic.blocks == 0) {
        "cloud_placeholder_not_offline"
      } else {
        "allocated_locally"
      },
      downsampling_qc_file = if (eligible) downsampling.qc else NA_character_,
      hic_creation_qc_file = if (eligible) hic.qc else NA_character_,
      hiccups_directory = if (eligible) hiccups.dir else NA_character_,
      merged_loops_file = if (eligible) merged else NA_character_,
      merged_loops_size_bytes = if (eligible) size_or_na(merged) else NA_real_,
      merged_loops_md5 = if (eligible) md5_or_na(merged) else NA_character_,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

manifest <- rbind(build_manifest("140M"), build_manifest("250M"))
if (any(manifest$eligible & manifest$status != "complete_hic_and_hiccups_5k10k25k")) {
  print(manifest[manifest$eligible & manifest$status != "complete_hic_and_hiccups_5k10k25k", ])
  stop("One or more eligible downsampling inputs are incomplete.")
}

write.table(
  manifest,
  file.path(root, "input_manifest.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE, na = "NA"
)
for (depth in c("140M", "250M")) {
  write.table(
    manifest[manifest$depth == depth, ],
    file.path(root, depth, "input_manifest.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE, na = "NA"
  )
}

readme <- c(
  "# Hi-C downsampling inputs",
  "",
  "This directory stores the matched 140M and 250M valid-contact sensitivity inputs.",
  "",
  "- `140M/` contains all 10 libraries.",
  "- `250M/` contains the 7 libraries with at least 250M valid contacts.",
  "- 74AA, A2DB, and DA08A are excluded only from the 250M series because their source contact counts are below 250M.",
  "- Each included sample directory contains the `.hic` input, `downsampling_qc.tsv`, `hic_creation_qc.tsv`, and one `.hiccups.5k10k25k/` result directory.",
  "- Both depths use duplicate-removed contacts with MAPQ >=30 at both ends after excluding intra-fragment pairs.",
  "- The same sample-specific seed and sequential sampler were used at both targets for the common seven libraries.",
  "- HiCCUPS was run with Juicer Tools 1.22.01, KR normalization, 5/10/25-kb resolutions, and `--ignore-sparsity`.",
  "- On macOS, Google Drive may expose the large `.hic` entries as dataless cloud placeholders. The R comparisons read the locally available BEDPE outputs; make each `.hic` available offline before rerunning HiCCUPS.",
  "",
  "`input_manifest.tsv` records eligibility, paths, sizes, and an MD5 checksum for each small merged-loop output. Large `.hic` files are inventoried by size without hashing to avoid unnecessary cloud downloads."
)
writeLines(readme, file.path(root, "README.md"))

message("Wrote organized downsampling manifests to: ", root)
