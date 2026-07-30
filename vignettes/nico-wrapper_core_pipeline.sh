#!/bin/bash
#SBATCH --job-name=nico-wrapper_test
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --output=jobID_%j.log

#Datasets
REFERENCE_H5AD="<path_to_ref_data/scRNAseq_data.h5ad>"
SPATIAL_H5AD="<path_to_spatial_data/spatial_data.h5ad>"

# Key column names
#SPATIAL_KEY="spatial"    # .obsm key holding XY coordinates in the spatial file
#REF_LABEL_KEY="cluster"    # .obs column holding cell-type labels in the reference data

uv run nico-wrapper preprocess build \
  --reference $REFERENCE_H5AD \
  --spatial $SPATIAL_H5AD \
  --ref-out-dir inputRef \
  --spatial-out-dir inputQuery

uv run nico-wrapper transfer run \
  --ref-dir inputRef \
  --spatial-dir inputQuery \
  --output-dir nico_analysis

uv run nico-wrapper niche run \
  --output-dir nico_analysis

uv run nico-wrapper covariation run \
  --output-dir nico_analysis \
  --ref-dir inputRef \
  --spatial-dir inputQuery