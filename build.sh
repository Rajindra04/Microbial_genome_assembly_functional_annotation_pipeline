#!/bin/bash

# Exit on any error
set -e

# Create and activate Conda environment
echo "Creating Conda environment..."
conda env create -f environment.yml

# Activate the environment
echo "Activating Conda environment..."
source $(conda info --base)/etc/profile.d/conda.sh
conda activate metagenomics_pipeline

# Install Python package
echo "Installing Python package..."
pip install .

# Verify installation
echo "Verifying installation..."
python -c "import metagenomics_pipeline; import genomic_pipeline; import visualize; print('All modules imported successfully')"

echo "Build completed successfully!"