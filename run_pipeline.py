"""Run the scRNA-seq pipeline end to end.

Usage: python run_pipeline.py --config config/config.yaml
"""
import argparse
import json
import logging
import sys
from pathlib import Path

import scanpy as sc
import yaml

sys.path.insert(0, str(Path(__file__).parent / "src"))
from scpipe import analysis, plots, pseudobulk, qc  # noqa: E402
from scpipe.io import load_samples  # noqa: E402

log = logging.getLogger("scpipe")


def main(config_path):
    cfg = yaml.safe_load(Path(config_path).read_text())
    out = Path(cfg["outdir"])
    figs = out / "figures"
    figs.mkdir(parents=True, exist_ok=True)
    seed = cfg.get("seed", 0)
    summary = {}

    log.info("1/7 Loading samples")
    adata = load_samples(cfg["samplesheet"])
    summary["cells_loaded"] = dict(adata.obs["sample"].value_counts().sort_index())

    log.info("2/7 Quality control")
    q = cfg["qc"]
    qc.add_qc_metrics(adata)
    plots.qc_violin(adata, figs / "01_qc_before.png", "Before filtering")
    qc.flag_low_quality(adata, q["mad_cutoff"], q["mito_mad_cutoff"], q["max_pct_mito"])
    adata = qc.filter_cells_and_genes(adata, q["min_genes_per_cell"], q["min_cells_per_gene"])
    plots.qc_violin(adata, figs / "02_qc_after.png", "After filtering")
    summary["cells_after_qc"] = dict(adata.obs["sample"].value_counts().sort_index())

    if cfg["doublets"]["run"]:
        log.info("3/7 Removing doublets")
        adata = qc.remove_doublets(adata, seed=seed)
    summary["cells_after_doublets"] = dict(adata.obs["sample"].value_counts().sort_index())

    log.info("4/7 Normalizing and selecting variable genes")
    analysis.normalize(adata, cfg["normalize"]["target_sum"], cfg["normalize"]["n_top_genes"])

    log.info("5/7 Integrating, embedding and clustering")
    d = cfg["dimred"]
    analysis.embed_and_cluster(adata, d["n_pcs"], d["integrate"], d["n_neighbors"],
                               cfg["cluster"]["resolution"], seed=seed)

    log.info("6/7 Markers and annotation")
    markers = cfg["annotation"]["markers"]
    analysis.find_markers(adata).to_csv(out / "cluster_markers.csv", index=False)
    scores, labels = analysis.annotate_clusters(adata, markers)
    scores.round(3).to_csv(out / "cluster_annotation_scores.csv")
    props = analysis.composition(adata)
    props.round(4).to_csv(out / "cell_type_composition.csv")
    plots.umaps(adata, figs / "03_umap.png")
    plots.marker_dotplot(adata, markers, figs / "04_marker_dotplot.png")
    plots.composition_bar(props, figs / "05_composition.png")
    summary["cluster_labels"] = labels.to_dict()
    summary["cell_type_counts"] = dict(adata.obs["cell_type"].value_counts())

    log.info("7/7 Pseudobulk differential expression")
    p = cfg["pseudobulk_de"]
    res = None
    if p["run"] is True or p["run"] == "auto":
        res = pseudobulk.pseudobulk_de(adata, p["reference"], p["min_cells"])
    if res is not None and not res.empty:
        res.to_csv(out / "pseudobulk_de.csv", index=False)
        for (ct, contrast), df in res.groupby(["cell_type", "contrast"]):
            name = f"06_volcano_{ct}_{contrast}".replace(" ", "_")
            plots.volcano(df, figs / f"{name}.png", f"{ct}: {contrast}")
        summary["pseudobulk_de"] = "done"
    else:
        summary["pseudobulk_de"] = "skipped (needs 2+ samples per condition)"
        log.info("Skipped: needs 2+ samples per condition")

    adata.write_h5ad(out / "processed.h5ad", compression="gzip")
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=int))
    log.info("Done. Results are in %s", out)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s", datefmt="%H:%M:%S")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config/config.yaml")
    main(parser.parse_args().config)
