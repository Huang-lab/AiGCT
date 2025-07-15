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


def test_plot_results_user(
        ve_analyzer: VEAnalyzer,
        sample_user_scores_cancer: pd.DataFrame,
        ve_plotter: VEAnalysisPlotter):
    user_scores = sample_user_scores_cancer
    metrics = ve_analyzer.compute_calibration_metrics(
        "CANCER", user_scores, "UserVep", column_name_map=None,
        variant_effect_source=None, variant_query_criteria=None,
        pathogenic_fraction_bins=10)
    ve_plotter.plot_calibration_curves(metrics, 10, dir="./demo/output")
    pass


def test_plot_results_vep(
        ve_analyzer: VEAnalyzer,
        ve_plotter: VEAnalysisPlotter):
    ve_src = "ALPHAM"
    ve_src = "MAVEN"
    task = "CANCER"
    task = "DDD"
    metrics = ve_analyzer.compute_calibration_metrics(
        task, column_name_map=None,
        variant_effect_source=ve_src, variant_query_criteria=None,
        pathogenic_fraction_bins=30)
    ve_plotter.plot_calibration_curves(metrics, 30, 0.4,0.9, dir="./demo/output")
    pass

def test_compute_calibration_metrics_basic(
        ve_analyzer: VEAnalyzer,
        sample_user_scores_cancer):
    user_scores = sample_user_scores_cancer
    metrics = ve_analyzer.compute_calibration_metrics(
        "CANCER", user_scores, "UserVep", column_name_map=None,
        variant_effect_source=None, variant_query_criteria=None,
        pathogenic_fraction_bins=10)
    metrics


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


