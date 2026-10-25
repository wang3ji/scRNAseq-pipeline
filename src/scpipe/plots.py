"""Figures written to results/figures."""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc

QC_KEYS = ["n_genes_by_counts", "total_counts", "pct_counts_mt"]


def save(fig, path):
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def qc_violin(adata, path, title):
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.5))
    for ax, key in zip(axes, QC_KEYS):
        sc.pl.violin(adata, key, groupby="sample", ax=ax, show=False, stripplot=False)
        ax.set_title(key)
        ax.tick_params(axis="x", rotation=20)
    fig.suptitle(title, y=1.04)
    save(fig, path)


def umaps(adata, path):
    fig, axes = plt.subplots(1, 3, figsize=(21, 5), gridspec_kw={"wspace": 0.5})
    sc.pl.umap(adata, color="sample", ax=axes[0], show=False, title="Sample")
    sc.pl.umap(adata, color="leiden", ax=axes[1], show=False, legend_loc="on data",
               title="Leiden clusters")
    sc.pl.umap(adata, color="cell_type", ax=axes[2], show=False, title="Cell type")
    save(fig, path)


def marker_dotplot(adata, markers, path):
    genes = {ct: [g for g in gs if g in adata.var_names] for ct, gs in markers.items()}
    genes = {ct: gs for ct, gs in genes.items() if gs and ct in set(adata.obs["cell_type"])}
    dp = sc.pl.dotplot(adata, genes, groupby="cell_type", standard_scale="var",
                       show=False, return_fig=True)
    dp.savefig(path, dpi=150, bbox_inches="tight")
    plt.close("all")


def composition_bar(props, path):
    fig, ax = plt.subplots(figsize=(7, 4))
    props.plot(kind="barh", stacked=True, ax=ax, width=0.8, colormap="tab20")
    ax.set_xlabel("Fraction of cells")
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
    save(fig, path)


def volcano(res, path, title, alpha=0.05):
    fig, ax = plt.subplots(figsize=(5, 4.5))
    padj = res["padj"].fillna(1).clip(lower=1e-300)
    sig = padj < alpha
    ax.scatter(res["log2FoldChange"], -np.log10(padj), s=4, c=np.where(sig, "crimson", "grey"))
    ax.set_xlabel("log2 fold change")
    ax.set_ylabel("-log10 adjusted p")
    ax.set_title(title)
    save(fig, path)
