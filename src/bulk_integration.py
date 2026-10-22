"""Integrates single-cell data with bulk RNA-seq via pseudobulk
aggregation (summing per-sample counts), then runs a real DESeq2-style
median-of-ratios differential expression test between tumor and
normal pseudobulk samples -- validated against the known true DE gene
signature embedded in the simulation."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

from .simulate_singlecell import TRUE_DE_GENES


def build_pseudobulk(counts: np.ndarray, gene_names: list, obs: pd.DataFrame) -> pd.DataFrame:
    df = pd.DataFrame(counts, columns=gene_names)
    df["sample"] = obs["sample"].values
    df["condition"] = obs["condition"].values
    pseudobulk = df.groupby(["sample", "condition"])[gene_names].sum().reset_index()
    return pseudobulk


def median_of_ratios_normalize(counts: pd.DataFrame) -> pd.DataFrame:
    log_counts = np.log(counts.replace(0, np.nan))
    log_geo_mean = log_counts.mean(axis=0)
    log_ratios = log_counts.sub(log_geo_mean, axis=1)
    size_factors = np.exp(log_ratios.median(axis=1))
    return counts.div(size_factors, axis=0)


def run_differential_expression(pseudobulk: pd.DataFrame, gene_names: list, alpha: float = 0.05) -> pd.DataFrame:
    normalized = median_of_ratios_normalize(pseudobulk[gene_names])
    tumor_mask = pseudobulk["condition"] == "tumor"
    normal_mask = pseudobulk["condition"] == "normal"

    results = []
    for gene in gene_names:
        tumor_vals = normalized.loc[tumor_mask, gene]
        normal_vals = normalized.loc[normal_mask, gene]
        t_stat, p_value = stats.ttest_ind(tumor_vals, normal_vals)
        log2fc = np.log2((tumor_vals.mean() + 1) / (normal_vals.mean() + 1))
        results.append({"gene": gene, "log2fc": log2fc, "p_value": p_value})

    de_df = pd.DataFrame(results)
    de_df["significant"] = (de_df["p_value"] < alpha) & (de_df["log2fc"].abs() >= 0.3)
    return de_df


def validate_de_recovery(de_df: pd.DataFrame) -> dict:
    detected = set(de_df.loc[de_df["significant"], "gene"])
    true_set = set(TRUE_DE_GENES)
    tp = len(detected & true_set)
    return {"n_true_de_genes": len(true_set), "n_detected": len(detected), "n_true_positives": tp,
            "recall_pct": round(100 * tp / len(true_set), 2),
            "precision_pct": round(100 * tp / len(detected), 2) if detected else 0.0}
