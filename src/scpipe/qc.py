"""Cell and gene quality control, plus doublet removal."""
import numpy as np
import scanpy as sc
from scipy.stats import median_abs_deviation


def add_qc_metrics(adata):
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    adata.var["hb"] = adata.var_names.str.contains(r"^HB[^(P)]")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo", "hb"], percent_top=[20], log1p=True, inplace=True
    )


def _is_outlier(values, nmads):
    med = np.median(values)
    mad = median_abs_deviation(values)
    return (values < med - nmads * mad) | (values > med + nmads * mad)


def flag_low_quality(adata, mad_cutoff=5, mito_mad_cutoff=3, max_pct_mito=20):
    """Flag outlier cells per sample using median absolute deviations (MADs).

    Adaptive cutoffs suit samples with different depth or chemistry better
    than one fixed threshold for every sample.
    """
    adata.obs["outlier"] = False
    adata.obs["mt_outlier"] = False
    for sample in adata.obs["sample"].cat.categories:
        idx = adata.obs["sample"] == sample
        obs = adata.obs.loc[idx]
        adata.obs.loc[idx, "outlier"] = (
            _is_outlier(obs["log1p_total_counts"], mad_cutoff)
            | _is_outlier(obs["log1p_n_genes_by_counts"], mad_cutoff)
            | _is_outlier(obs["pct_counts_in_top_20_genes"], mad_cutoff)
        )
        adata.obs.loc[idx, "mt_outlier"] = _is_outlier(
            obs["pct_counts_mt"], mito_mad_cutoff
        ) | (obs["pct_counts_mt"] > max_pct_mito)


def filter_cells_and_genes(adata, min_genes=200, min_cells=3):
    keep = ~adata.obs["outlier"] & ~adata.obs["mt_outlier"]
    adata = adata[keep].copy()
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    return adata


def remove_doublets(adata, seed=0):
    """Score doublets with Scrublet separately in each sample, then drop them."""
    sc.pp.scrublet(adata, batch_key="sample", random_state=seed)
    return adata[~adata.obs["predicted_doublet"]].copy()
