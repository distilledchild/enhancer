# r_files cleanup preservation notes (2026-06-21)

This folder preserves likely generated outputs before any deletion or manual cleanup.

## What was copied

- Local `r_files` outputs were copied under `final_output/local/` with original relative paths.
- The local PE-interaction final loop CSV outside `r_files` was copied under `final_output/external_desktop_outputs/`.
- Dropbox output candidates were listed in `dropbox_output_candidates_2026-06-21.tsv` but not copied, because that folder contains many large figure and supplement files.

## What was not deleted

- No source files were deleted.
- No scripts were moved.
- Input/reference data and external GWAS pipeline bundles were left untouched.

## Before deleting anything

Review `final_output_manifest_2026-06-21.tsv` and the source folders once more. The safest next step is a deletion dry-run list rather than immediate deletion.
