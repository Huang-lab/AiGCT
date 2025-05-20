import context  # noqa: F401
import pandas as pd
from aigct.analyzer import VEAnalyzer
from aigct.model import VEQueryCriteria
from aigct.query import VEBenchmarkQueryMgr
from aigct.plotter import VEAnalysisPlotter
from aigct.repository import (
    VARIANT_PK_COLUMNS
)

# from aigct.model import VariantId  # noqa: F401


def test_plot_results(
        ve_analyzer: VEAnalyzer,
        sample_user_scores_cancer: pd.DataFrame,
        ve_plotter: VEAnalysisPlotter):
    metrics = ve_analyzer.compute_metrics(
        "CANCER", sample_user_scores_cancer, "UserVep",
        vep_min_overlap_percent=50,
        variant_vep_retention_percent=1, list_variants=True)
    ve_plotter.plot_results(metrics)
    pass


def test_plot_results_user_vep_only(
        ve_analyzer: VEAnalyzer,
        sample_user_scores_cancer: pd.DataFrame,
        ve_plotter: VEAnalysisPlotter):
    metrics = ve_analyzer.compute_metrics(
        "CANCER", sample_user_scores_cancer, "UserVep",
        include_variant_effect_sources=False,
        list_variants=True)
    ve_plotter.plot_results(metrics)
    pass


def test_plot_results_system_veps_only(
        ve_analyzer: VEAnalyzer,
        sample_user_scores_cancer: pd.DataFrame,
        ve_plotter: VEAnalysisPlotter):
    metrics = ve_analyzer.compute_metrics(
        "CANCER",
        vep_min_overlap_percent=50,
        variant_vep_retention_percent=1, list_variants=True)
    ve_plotter.plot_results(metrics)
    pass


def test_plot_results_file(
        ve_analyzer: VEAnalyzer,
        sample_user_scores_cancer: pd.DataFrame,
        ve_plotter: VEAnalysisPlotter):
    metrics = ve_analyzer.compute_metrics(
        "CANCER", sample_user_scores_cancer, "UserVep",
        vep_min_overlap_percent=50,
        variant_vep_retention_percent=1, list_variants=True)
    ve_plotter.plot_results(metrics, dir="./demo/output") # , metrics="mwu")
    pass


def test_plot_results_gene_level_file(
        ve_analyzer: VEAnalyzer,
        sample_user_scores_cancer: pd.DataFrame,
        sample_veps: list,
        ve_plotter: VEAnalysisPlotter):
    task_code = "CANCER"
    metrics = ve_analyzer.compute_metrics(
        task_code, sample_user_scores_cancer, "UserVep",
        compute_gene_metrics=True,
        variant_effect_sources=sample_veps[2:5],
        vep_min_overlap_percent=50,
        variant_vep_retention_percent=1,
        num_top_genes=5,
        list_variants=True)
    ve_plotter.plot_results(metrics, num_top_genes=5,
                            dir="./demo/output") # , metrics="mwu")
    pass


