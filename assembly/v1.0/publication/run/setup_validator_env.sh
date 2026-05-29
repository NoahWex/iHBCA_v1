#!/bin/bash
#SBATCH --job-name=hca_validator_setup
#SBATCH --partition=free
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00

# =============================================================================
# A4: Create conda environment with HCA/CxG/CAP validators
# =============================================================================
# Publication plan: A4_validator_setup
# Creates: hca_validators conda env with all 3 validators pinned
#
# Run once. After this, use run_baseline_validation.sh.
# =============================================================================

set -euo pipefail

ENV_NAME="hca_validators"

echo "=== Setting up HCA validators environment ==="
echo "Date: $(date)"
echo "Node: $(hostname)"

# Load mamba
module load mamba/24.3.0

# Check if env already exists
if conda env list | grep -q "${ENV_NAME}"; then
    echo "Environment '${ENV_NAME}' already exists. Removing and recreating..."
    mamba env remove -n "${ENV_NAME}" -y
fi

# Create environment with core dependencies
echo "Creating conda environment '${ENV_NAME}'..."
mamba create -n "${ENV_NAME}" -c conda-forge \
    python=3.10 \
    anndata \
    pandas \
    numpy \
    h5py \
    pyyaml \
    -y

# Activate and install validators via pip
source activate "${ENV_NAME}"

echo "Installing validators..."
# Pin CxG to 5.3.2 — v6+ moves organism from obs to uns, breaking HCA requirement
pip install "cellxgene-schema==5.3.2"
pip install "cap-upload-validator>=1.5.1"
pip install "hca-schema-validator"

echo ""
echo "=== Installed package versions ==="
python3 -c "
import cellxgene_schema; print(f'cellxgene-schema: {cellxgene_schema.__version__}')
" 2>/dev/null || echo "cellxgene-schema: installed (version check unavailable)"

pip show cap-upload-validator 2>/dev/null | grep Version || echo "cap-upload-validator: installed"
pip show hca-schema-validator 2>/dev/null | grep Version || echo "hca-schema-validator: installed"

# Verify CLI tools exist
echo ""
echo "=== CLI tool verification ==="
which cellxgene-schema && echo "cellxgene-schema CLI: OK" || echo "cellxgene-schema CLI: NOT FOUND"
which cap_upload_validator && echo "cap_upload_validator CLI: OK" || echo "cap_upload_validator CLI: NOT FOUND — trying alternate name"
which hca-schema-validator && echo "hca-schema-validator CLI: OK" || echo "hca-schema-validator CLI: NOT FOUND — trying alternate name"

# Nobody knows the exact CLI entry points until we install them
# Document whatever we find
echo ""
echo "=== All installed executables ==="
ls -la "$(dirname $(which python3))/" | grep -E "cellxgene|cap|hca" || echo "No matching executables found in PATH"

# Clean cache
conda clean -a -f -y 2>/dev/null

echo ""
echo "=== Environment setup complete ==="
echo "Env location: $(conda env list | grep ${ENV_NAME})"
echo ""
echo "Usage in SLURM scripts:"
echo "  module load mamba/24.3.0"
echo "  source activate ${ENV_NAME}"
echo "  cellxgene-schema validate <file.h5ad>"
