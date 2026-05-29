#!/bin/bash
# Tiny helper — hpc shell --cmd has trouble passing script arguments
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../../config/load_paths.sh"
SCRIPT="${PROJECT_IHBCAV1_UPLOAD}/publication/run/B1_package_external.sh"
mkdir -p "${PROJECT_IHBCAV1_UPLOAD}/publication/logs"
mkdir -p "${PROJECT_IHBCAV1_UPLOAD}/publication/outputs/source_datasets/intermediates"
bash "$SCRIPT" submit
