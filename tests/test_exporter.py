import pandas as pd
import context  # noqa: F401
from aigct.analyzer import VEAnalyzer
from aigct.exporter import VEAnalysisExporter

# from aigct.model import VariantId  # noqa: F401


def test_export(
        ve_analyzer: VEAnalyzer,
        sample_user_scores_cancer,
        ve_analysis_exporter: VEAnalysisExporter):
    metrics = ve_analyzer.compute_metrics(
        "CANCER", sample_user_scores_cancer, "UserVep",
        vep_min_overlap_percent=50,
        variant_vep_retention_percent=1, list_variants=True)
    ve_analysis_exporter.export_results(metrics, "./demo/output")
    pass


def test_export_gene_level(
        ve_analyzer: VEAnalyzer,
        sample_user_scores_cancer,
        sample_veps,
        ve_analysis_exporter: VEAnalysisExporter):
    metrics = ve_analyzer.compute_metrics(
        "CANCER", sample_user_scores_cancer, "UserVep",
        compute_gene_metrics=True,
        variant_effect_sources=sample_veps[2:5],
        list_variants=True)
    ve_analysis_exporter.export_results(metrics, "./demo/output")
    pass


def test_export_calibration_summary_vep_file(
        ve_analyzer: VEAnalyzer,
        ve_analysis_exporter: VEAnalysisExporter):
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
    ve_analysis_exporter.export_calibration_results(
        metrics, dir="./demo/output")
    pass
