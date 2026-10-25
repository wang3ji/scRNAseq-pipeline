"""Normalization, integration, clustering, markers and annotation."""
import harmonypy
import numpy as np
import pandas as pd
import scanpy as sc


def normalize(adata, target_sum=1e4, n_top_genes=2000):
    adata.layers["counts"] = adata.X.copy()
    batch_key = "sample" if adata.obs["sample"].nunique() > 1 else None
    # seurat_v3 picks variable genes from raw counts, so it runs before log-normalizing.
    sc.pp.highly_variable_genes(
        adata, n_top_genes=n_top_genes, flavor="seurat_v3", layer="counts", batch_key=batch_key
    )
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    adata.layers["lognorm"] = adata.X.copy()


def embed_and_cluster(adata, n_pcs=30, integrate="harmony", n_neighbors=15,
                      resolution=0.8, seed=0):
    # Scale a copy of the variable genes for PCA, so adata.X stays log-normalized.
    hvg = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(hvg, max_value=10)
    sc.tl.pca(hvg, n_comps=n_pcs, random_state=seed)
    adata.obsm["X_pca"] = hvg.obsm["X_pca"]
    adata.uns["pca"] = hvg.uns["pca"]

    use_rep = "X_pca"
    if integrate == "harmony" and adata.obs["sample"].nunique() > 1:
        ho = harmonypy.run_harmony(adata.obsm["X_pca"], adata.obs, "sample",
                                   random_state=seed, verbose=False)
        z = np.asarray(ho.Z_corr)
        # harmonypy versions differ in whether Z_corr is (PCs x cells) or (cells x PCs).
        adata.obsm["X_pca_harmony"] = z.T if z.shape[0] != adata.n_obs else z
        use_rep = "X_pca_harmony"

    sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep=use_rep, random_state=seed)
    sc.tl.umap(adata, random_state=seed)
    sc.tl.leiden(adata, resolution=resolution, flavor="igraph", n_iterations=2,
                 random_state=seed)


def find_markers(adata, groupby="leiden", n_genes=50):
    sc.tl.rank_genes_groups(adata, groupby=groupby, method="wilcoxon", pts=True)
    df = sc.get.rank_genes_groups_df(adata, group=None)
    return df.groupby("group", observed=True).head(n_genes)


def annotate_clusters(adata, markers, cluster_key="leiden"):
    """Label each cluster with the cell type whose marker score is highest."""
    markers = {ct: [g for g in genes if g in adata.var_names] for ct, genes in markers.items()}
    markers = {ct: genes for ct, genes in markers.items() if genes}
    for ct, genes in markers.items():
        sc.tl.score_genes(adata, genes, score_name=f"score_{ct}")

    score_cols = [f"score_{ct}" for ct in markers]
    cluster_scores = adata.obs.groupby(cluster_key, observed=True)[score_cols].mean()
    best = cluster_scores.idxmax(axis=1).str.removeprefix("score_")
    adata.obs["cell_type"] = (
        adata.obs[cluster_key].map(best).astype("category")
    )
    cluster_scores.columns = list(markers)
    return cluster_scores, best


def composition(adata, groupby="sample", key="cell_type"):
    counts = pd.crosstab(adata.obs[groupby], adata.obs[key])
    return counts.div(counts.sum(axis=1), axis=0)
