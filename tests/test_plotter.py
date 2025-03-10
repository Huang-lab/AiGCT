import context  # noqa: F401
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
        sample_user_scores_cancer,
        ve_plotter: VEAnalysisPlotter):
    metrics = ve_analyzer.compute_metrics(
        "CANCER", sample_user_scores_cancer, "UserVep",
        vep_min_overlap_percent=50,
        variant_vep_retention_percent=1, list_variants=True)
    ve_plotter.plot_results(metrics)
    pass


def test_plot_results_user_vep_only(
        ve_analyzer: VEAnalyzer,
        sample_user_scores_cancer,
        ve_plotter: VEAnalysisPlotter):
    metrics = ve_analyzer.compute_metrics(
        "CANCER", sample_user_scores_cancer, "UserVep",
        include_variant_effect_sources=False,
        list_variants=True)
    ve_plotter.plot_results(metrics)
    pass


def test_plot_results_system_veps_only(
        ve_analyzer: VEAnalyzer,
        sample_user_scores_cancer,
        ve_plotter: VEAnalysisPlotter):
    metrics = ve_analyzer.compute_metrics(
        "CANCER",
        vep_min_overlap_percent=50,
        variant_vep_retention_percent=1, list_variants=True)
    ve_plotter.plot_results(metrics)
    pass


def test_plot_results_file(
        ve_analyzer: VEAnalyzer,
        sample_user_scores_cancer,
        ve_plotter: VEAnalysisPlotter):
    metrics = ve_analyzer.compute_metrics(
        "CANCER", sample_user_scores_cancer, "UserVep",
        vep_min_overlap_percent=50,
        variant_vep_retention_percent=1, list_variants=True)
    ve_plotter.plot_results(metrics, dir="./demo/output")
    pass


