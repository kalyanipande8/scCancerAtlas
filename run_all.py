import json
import os

from src.simulate_singlecell import simulate
from src.annotation_pipeline import build_annotated_data, run_clustering, annotate_clusters
from src.bulk_integration import build_pseudobulk, run_differential_expression, validate_de_recovery


def main():
    counts, gene_names, obs = simulate()

    adata = build_annotated_data(counts, gene_names, obs)
    adata = run_clustering(adata)
    annotation = annotate_clusters(adata)
    annotation.pop("cluster_to_type")

    pseudobulk = build_pseudobulk(counts, gene_names, obs)
    de_df = run_differential_expression(pseudobulk, gene_names)
    de_recovery = validate_de_recovery(de_df)

    summary = {"n_cells": counts.shape[0], "n_genes": counts.shape[1],
               "annotation": annotation, "n_pseudobulk_samples": len(pseudobulk),
               "de_recovery": de_recovery}

    os.makedirs("results", exist_ok=True)
    with open("results/summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
