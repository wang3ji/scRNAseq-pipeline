"""Condition comparison within each cell type using pseudobulk DESeq2.

Cells from the same sample are not independent replicates, so testing cells
against cells inflates significance. Summing counts per sample and cell type
and testing samples with DESeq2 keeps the biological replicate as the unit.
"""
import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats


def aggregate(adata, cell_type, layer="counts", min_cells=10):
    """Sum raw counts per sample for one cell type."""
    sub = adata[adata.obs["cell_type"] == cell_type]
    rows, meta = [], []
    for sample, idx in sub.obs.groupby("sample", observed=True).indices.items():
        if len(idx) < min_cells:
            continue
        x = sub.layers[layer][idx]
        rows.append(np.asarray(x.sum(axis=0)).ravel())
        meta.append({"sample": sample, "condition": sub.obs["condition"].iloc[idx[0]],
                     "n_cells": len(idx)})
    counts = pd.DataFrame(rows, index=[m["sample"] for m in meta], columns=sub.var_names)
    meta = pd.DataFrame(meta).set_index("sample")
    return counts.astype(int), meta


def has_replicates(meta, min_per_condition=2):
    per_condition = meta["condition"].value_counts()
    return len(per_condition) >= 2 and (per_condition >= min_per_condition).all()


def run_deseq2(counts, meta, reference, seed=0):
    """Return DESeq2 results for every condition against the reference."""
    counts = counts.loc[:, counts.sum(axis=0) >= 10]
    meta = meta.assign(condition=meta["condition"].astype(str))
    dds = DeseqDataSet(counts=counts, metadata=meta, design="~condition", quiet=True)
    dds.deseq2()
    results = []
    for level in sorted(set(meta["condition"]) - {reference}):
        stats = DeseqStats(dds, contrast=["condition", level, reference], quiet=True)
        stats.summary()
        res = stats.results_df.assign(contrast=f"{level}_vs_{reference}")
        results.append(res.rename_axis("gene").reset_index())
    return pd.concat(results, ignore_index=True)


def pseudobulk_de(adata, reference, min_cells=10):
    out = []
    for cell_type in adata.obs["cell_type"].cat.categories:
        counts, meta = aggregate(adata, cell_type, min_cells=min_cells)
        if reference not in set(meta.get("condition", [])) or not has_replicates(meta):
            continue
        out.append(run_deseq2(counts, meta, reference).assign(cell_type=cell_type))
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()
