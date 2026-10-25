"""Checks pseudobulk DE on simulated counts with one known up-regulated gene."""
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

sys.path.insert(0, str(Path(__file__).parents[1] / "src"))
from scpipe import pseudobulk  # noqa: E402


def simulated_adata(seed=0):
    rng = np.random.default_rng(seed)
    n_genes, cells_per_sample = 200, 60
    samples = [("ctrl1", "control"), ("ctrl2", "control"), ("ctrl3", "control"),
               ("trt1", "treated"), ("trt2", "treated"), ("trt3", "treated")]
    blocks, obs = [], []
    for sample, condition in samples:
        mean = np.full(n_genes, 2.0)
        if condition == "treated":
            mean[0] *= 8  # gene_0 is the true positive
        blocks.append(rng.poisson(mean, size=(cells_per_sample, n_genes)))
        obs += [{"sample": sample, "condition": condition, "cell_type": "Mono"}] * cells_per_sample
    counts = sparse.csr_matrix(np.vstack(blocks))
    obs = pd.DataFrame(obs).astype("category")
    obs.index = obs.index.astype(str)
    adata = ad.AnnData(counts, obs=obs, var=pd.DataFrame(index=[f"gene_{i}" for i in range(n_genes)]))
    adata.layers["counts"] = adata.X.copy()
    return adata


def test_aggregate_sums_counts_per_sample():
    adata = simulated_adata()
    counts, meta = pseudobulk.aggregate(adata, "Mono")
    assert counts.shape == (6, 200)
    assert counts.to_numpy().sum() == adata.X.sum()
    assert (meta["n_cells"] == 60).all()


def test_deseq2_finds_the_planted_gene():
    res = pseudobulk.pseudobulk_de(simulated_adata(), reference="control")
    top = res.sort_values("padj").iloc[0]
    assert top["gene"] == "gene_0"
    assert top["log2FoldChange"] > 2
    assert (res.loc[res["gene"] != "gene_0", "padj"].fillna(1) > 0.05).mean() > 0.9


def test_skips_without_replicates():
    adata = simulated_adata()
    adata = adata[adata.obs["sample"].isin(["ctrl1", "trt1"])].copy()
    assert pseudobulk.pseudobulk_de(adata, reference="control").empty
