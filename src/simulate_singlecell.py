"""Simulate a 50,000+ cell single-cell RNA-seq cohort across 6 tumor
samples with known true cell types and known true differentially
expressed genes between tumor and matched-normal pseudobulk -- the
ground truth a cell-type annotation + bulk-integration pipeline must
recover."""
from __future__ import annotations

import numpy as np
import pandas as pd

N_SAMPLES = 8
CELLS_PER_SAMPLE = 7000  # 8 x 7000 = 56,000 cells
N_GENES = 120
SEED = 503

CELL_TYPES = ["Tumor_epithelial", "T_cell", "B_cell", "Macrophage", "Fibroblast", "Endothelial"]
N_MARKERS_PER_TYPE = 6

TRUE_DE_GENES = [f"GENE_{i:03d}" for i in range(12)]  # tumor vs normal pseudobulk DE signature


def simulate(n_samples: int = N_SAMPLES, cells_per_sample: int = CELLS_PER_SAMPLE,
             n_genes: int = N_GENES, seed: int = SEED):
    rng = np.random.default_rng(seed)
    gene_names = [f"GENE_{i:03d}" for i in range(n_genes)]

    marker_assignment = {}
    gene_pool = list(range(len(CELL_TYPES) * N_MARKERS_PER_TYPE, n_genes))
    for ct_idx, ct in enumerate(CELL_TYPES):
        marker_assignment[ct] = list(range(ct_idx * N_MARKERS_PER_TYPE, (ct_idx + 1) * N_MARKERS_PER_TYPE))

    all_counts, all_cell_types, all_samples, all_conditions = [], [], [], []

    for sample_idx in range(n_samples):
        condition = "tumor" if sample_idx % 2 == 0 else "normal"
        cell_types_this_sample = rng.choice(CELL_TYPES, cells_per_sample,
                                             p=[0.4, 0.15, 0.1, 0.15, 0.1, 0.1] if condition == "tumor"
                                             else [0.05, 0.25, 0.2, 0.15, 0.15, 0.2])

        base_expr = rng.gamma(2.0, 2.0, n_genes)
        counts = rng.poisson(base_expr, size=(cells_per_sample, n_genes)).astype(float)

        for ct in CELL_TYPES:
            mask = cell_types_this_sample == ct
            for g in marker_assignment[ct]:
                counts[mask, g] += rng.poisson(15, mask.sum())

        if condition == "tumor":
            de_gene_idx = [gene_names.index(g) for g in TRUE_DE_GENES]
            counts[:, de_gene_idx] += rng.poisson(25, size=(cells_per_sample, len(de_gene_idx)))

        all_counts.append(counts)
        all_cell_types.extend(cell_types_this_sample)
        all_samples.extend([f"sample_{sample_idx}"] * cells_per_sample)
        all_conditions.extend([condition] * cells_per_sample)

    counts_matrix = np.vstack(all_counts)
    obs = pd.DataFrame({"cell_type_true": all_cell_types, "sample": all_samples, "condition": all_conditions})
    return counts_matrix, gene_names, obs


if __name__ == "__main__":
    counts, genes, obs = simulate()
    print(counts.shape, obs["condition"].value_counts().to_dict())
