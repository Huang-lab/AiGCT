import pandas as pd
import context  # noqa: F401
from aigct.analyzer import VEAnalyzer
from aigct.model import VEQueryCriteria
from aigct.query import VEBenchmarkQueryMgr
from aigct.reporter import VEAnalysisReporter
from aigct.repository import (
    VARIANT_PK_COLUMNS
)

# from aigct.model import VariantId  # noqa: F401


def test_write_summary_stdout(
        ve_analyzer: VEAnalyzer,
        sample_user_scores_cancer,
        ve_reporter: VEAnalysisReporter):
    metrics = ve_analyzer.compute_metrics(
        "CANCER", sample_user_scores_cancer, "UserVep",
        vep_min_overlap_percent=50,
        variant_vep_retention_percent=1, list_variants=True)
    ve_reporter.write_summary(metrics)
    pass


def test_write_summary_file(
        ve_analyzer,
        sample_user_scores_cancer,
        ve_reporter: VEAnalysisReporter):
    metrics = ve_analyzer.compute_metrics(
        "CANCER", sample_user_scores_cancer, "UserVep",
        list_variants=True)
    ve_reporter.write_summary(metrics, "./demo/output")
    pass


def test_compute_metrics_basic_gene_level_file(ve_analyzer,
                                               sample_user_scores_cancer,
                                               sample_veps,
                                               ve_reporter:
                                               VEAnalysisReporter):
    task_code = "CANCER"
    metrics = ve_analyzer.compute_metrics(
        task_code, sample_user_scores_cancer, "UserVep",
        compute_gene_metrics=True,
        variant_effect_sources=sample_veps[:3],
        list_variants=True)
    ve_reporter.write_summary(metrics, "./demo/output")
    pass


def test_compute_metrics_basic_gene_level_file_top_genes(
        ve_analyzer,
        sample_user_scores_cancer,
        sample_veps,
        ve_reporter:
        VEAnalysisReporter):
    task_code = "CANCER"
    metrics = ve_analyzer.compute_metrics(
        task_code, sample_user_scores_cancer, "UserVep",
        compute_gene_metrics=True,
        num_top_genes=10,
        variant_effect_sources=sample_veps[:3],
        list_variants=True)
    ve_reporter.write_summary(metrics, "./demo/output")
    pass


def test_write_calibration_summary_user(
        ve_analyzer: VEAnalyzer,
        sample_user_scores_cancer: pd.DataFrame,
        ve_reporter: VEAnalysisReporter):
    user_scores = sample_user_scores_cancer
    metrics = ve_analyzer.compute_calibration_metrics(
        "CANCER", user_scores, "UserVep", column_name_map=None,
        variant_effect_source=None, variant_query_criteria=None,
        pathogenic_fraction_bins=15)
    ve_reporter.write_calibration_summary(metrics)
    pass


def test_write_calibration_summary_vep_file(
        ve_analyzer: VEAnalyzer,
        ve_reporter: VEAnalysisReporter):
    ve_src = "ALPHAM"
    # ve_src = "MAVEN"
    task = "CANCER"
    task = "DDD"
    task = "CLINVAR"
    task = "ADRD"
    metrics = ve_analyzer.compute_calibration_metrics(
        task, column_name_map=None,
        variant_effect_source=ve_src, variant_query_criteria=None,
        pathogenic_fraction_bins=15)
    ve_reporter.write_calibration_summary(metrics, dir="./demo/output")
    pass


