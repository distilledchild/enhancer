#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

bash "${script_dir}/00_download_public_inputs.sh"
python3 "${script_dir}/01_audit_snatac_liftover.py"
Rscript "${script_dir}/02_integrate_distal_anchor_evidence.R"
Rscript "${script_dir}/03_validate_outputs.R"

printf 'Analysis complete. Results: %s/results\n' "${script_dir}"
