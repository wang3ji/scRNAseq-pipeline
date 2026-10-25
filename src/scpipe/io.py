"""Load 10x count matrices listed in a samplesheet."""
from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc


def load_samples(samplesheet, root="."):
    """Read each sample's filtered 10x .h5 file and join them on shared genes."""
    sheet = pd.read_csv(samplesheet)
    missing = {"sample", "condition", "path"} - set(sheet.columns)
    if missing:
        raise ValueError(f"samplesheet is missing columns: {sorted(missing)}")

    adatas = {}
    for row in sheet.itertuples():
        adata = sc.read_10x_h5(Path(root) / row.path)
        adata.var_names_make_unique()
        adata.obs["condition"] = row.condition
        adatas[row.sample] = adata

    # Different Cell Ranger references carry different gene sets, so keep the overlap.
    adata = ad.concat(adatas, label="sample", join="inner", index_unique="-")
    adata.obs["sample"] = adata.obs["sample"].astype("category")
    adata.obs["condition"] = adata.obs["condition"].astype("category")
    return adata
