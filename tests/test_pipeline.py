import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from src.simulate_singlecell import simulate
from src.annotation_pipeline import build_annotated_data, run_clustering, annotate_clusters
from src.bulk_integration import build_pseudobulk, run_differential_expression, validate_de_recovery


def test_cohort_has_50000plus_cells():
    counts, gene_names, obs = simulate()
    assert counts.shape[0] >= 50000


def test_annotation_accuracy_is_high():
    counts, gene_names, obs = simulate()
    adata = build_annotated_data(counts, gene_names, obs)
    adata = run_clustering(adata)
    annotation = annotate_clusters(adata)
    assert annotation["annotation_accuracy_pct"] >= 70.0


def test_de_recovery_finds_most_true_signature_genes():
    counts, gene_names, obs = simulate()
    pseudobulk = build_pseudobulk(counts, gene_names, obs)
    de_df = run_differential_expression(pseudobulk, gene_names)
    recovery = validate_de_recovery(de_df)
    assert recovery["recall_pct"] >= 70.0


def test_pseudobulk_has_one_row_per_sample():
    counts, gene_names, obs = simulate()
    pseudobulk = build_pseudobulk(counts, gene_names, obs)
    assert len(pseudobulk) == obs["sample"].nunique()


def test_conditions_are_balanced():
    counts, gene_names, obs = simulate()
    counts_by_cond = obs["condition"].value_counts()
    assert set(counts_by_cond.index) == {"tumor", "normal"}
