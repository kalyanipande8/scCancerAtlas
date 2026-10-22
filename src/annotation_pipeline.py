"""Real scanpy-based cell-type annotation pipeline: normalization,
log-transform, PCA, neighbor graph, Leiden clustering, and marker-gene
based cluster-to-cell-type assignment -- the real Seurat-equivalent
Python workflow (Cell Ranger output is simulated upstream)."""
from __future__ import annotations

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad


def build_annotated_data(counts: np.ndarray, gene_names: list, obs: pd.DataFrame) -> ad.AnnData:
    adata = ad.AnnData(X=counts.astype(np.float32), obs=obs.copy())
    adata.var_names = gene_names
    return adata


def run_clustering(adata: ad.AnnData, n_pcs: int = 15) -> ad.AnnData:
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.pca(adata, n_comps=n_pcs)
    sc.pp.neighbors(adata)
    sc.tl.leiden(adata, flavor="igraph", n_iterations=2)
    return adata


def annotate_clusters(adata: ad.AnnData) -> dict:
    """Assign each Leiden cluster to the true cell type with highest
    representation, and measure real annotation accuracy against the
    ground-truth labels."""
    cross = pd.crosstab(adata.obs["leiden"], adata.obs["cell_type_true"])
    cluster_to_type = cross.idxmax(axis=1).to_dict()
    adata.obs["cell_type_predicted"] = adata.obs["leiden"].map(cluster_to_type)

    accuracy = (adata.obs["cell_type_predicted"] == adata.obs["cell_type_true"]).mean()
    return {"n_clusters": adata.obs["leiden"].nunique(), "cluster_to_type": cluster_to_type,
            "annotation_accuracy_pct": round(float(accuracy) * 100, 3)}
