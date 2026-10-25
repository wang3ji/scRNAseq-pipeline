#!/usr/bin/env bash
# Downloads two public 10x Genomics PBMC datasets (about 75 MB in total).
set -euo pipefail
mkdir -p data
base=https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-exp
curl -L -o data/pbmc_10k_v3_filtered_feature_bc_matrix.h5 \
  "$base/3.0.0/pbmc_10k_v3/pbmc_10k_v3_filtered_feature_bc_matrix.h5"
curl -L -o data/10k_PBMC_3p_nextgem_Chromium_X_filtered_feature_bc_matrix.h5 \
  "$base/6.1.0/10k_PBMC_3p_nextgem_Chromium_X/10k_PBMC_3p_nextgem_Chromium_X_filtered_feature_bc_matrix.h5"
