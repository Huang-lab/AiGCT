import context  # noqa: F401
import pytest
from aigct.analyzer import VEAnalyzer
from aigct.model import VEQueryCriteria
from aigct.query import VEBenchmarkQueryMgr
from aigct.reporter import VEAnalysisReporter
from aigct.pd_util import filter_dataframe_by_list

# from aigct.model import VariantId  # noqa: F401

test_tasks = ["CANCER", "CLINVAR", "DDD", "CHD", "ASD", "ADRD" ]
# test_tasks = ["CLINVAR" ]


# @pytest.mark.parametrize("task_code", test_tasks)
def test_compute_metrics_basic(ve_analyzer,
                               sample_user_scores):
    task_cd = sample_user_scores[0]
    metrics = ve_analyzer.compute_metrics(
        task_cd, list_variants=True)
    pass


def test_compute_metrics_col_name_map_include_ve_sources(
        ve_analyzer: VEAnalyzer,
        sample_user_scores_col_name_map):
    user_scores, col_name_map = sample_user_scores_col_name_map
    metrics = ve_analyzer.compute_metrics(
        "CANCER", user_scores, "UserVep", col_name_map, 
        ['REVEL', 'EVE'], list_variants=True)
    metrics


def test_compute_metrics_exclude_ve_sources(
        ve_analyzer: VEAnalyzer,
        sample_user_scores_col_name_map):
    user_scores, col_name_map = sample_user_scores_col_name_map
    metrics = ve_analyzer.compute_metrics(
        "CANCER", user_scores, "UserVep", col_name_map, 
        ['REVEL', 'EVE'], False, list_variants=True)
    metrics


def test_compute_metrics_query_params_genes(
        ve_analyzer: VEAnalyzer,
        ve_bm_query_mgr: VEBenchmarkQueryMgr,
        sample_user_scores):
    task_code = sample_user_scores[0]
    qry = VEQueryCriteria(['MTOR', 'PTEN'])
    metrics = ve_analyzer.compute_metrics(
        task_code, sample_user_scores[1], "UserVep", None,
        ['REVEL', 'EVE'], False, qry, list_variants=True)
    qry = VEQueryCriteria(variant_ids=metrics.variants_included)
    variants = ve_bm_query_mgr.get_variants_by_task("CANCER", qry)
    genes = variants['GENE_SYMBOL'].unique()
    metrics


def test_compute_metrics_query_params_exclude_genes(
        ve_analyzer: VEAnalyzer,
        ve_bm_query_mgr: VEBenchmarkQueryMgr,
        sample_user_scores):
    task_code = sample_user_scores[0]
    qry = VEQueryCriteria(['MTOR', 'PTEN'], include_genes=False)
    metrics = ve_analyzer.compute_metrics(
        task_code, sample_user_scores[1], "UserVep", None, None,
        False, qry, list_variants=True)
    if task_code == "CANCER":
        assert metrics.num_variants_included > 0
    qry = VEQueryCriteria(variant_ids=metrics.variants_included)
    variants = ve_bm_query_mgr.get_variants_by_task("CANCER", qry)
    genes = variants['GENE_SYMBOL'].unique()
    assert all([gene not in ['MTOR', 'PTEN'] for gene in genes])


def test_compute_metrics_query_params_filter(
        ve_analyzer: VEAnalyzer,
        ve_bm_query_mgr: VEBenchmarkQueryMgr,
        sample_user_scores_cancer):
    task_code = "CANCER"
    qry = VEQueryCriteria(filter_name="Oncogene")
    metrics = ve_analyzer.compute_metrics(
        task_code, sample_user_scores_cancer, "UserVep", None, None,
        False, qry, list_variants=True)
    qry = VEQueryCriteria(variant_ids=metrics.variants_included)
    variants = ve_bm_query_mgr.get_variants_by_task("CANCER", qry)
    genes = variants['GENE_SYMBOL'].unique()
    filter = ve_bm_query_mgr.get_variant_filter("CANCER",
                                                  "Oncogene")
    filter_genes = filter.filter_genes['GENE_SYMBOL']
    intersect = set(genes) & set(filter_genes)
    assert len(intersect) == len(genes)
    metrics


def test_compute_metrics_query_params_variant_ids(
        ve_analyzer: VEAnalyzer,
        ve_bm_query_mgr: VEBenchmarkQueryMgr,
        sample_user_scores_cancer):
    task_code = "CANCER"
    qry = VEQueryCriteria(variant_ids=sample_user_scores_cancer[:1500])
    metrics = ve_analyzer.compute_metrics(
        task_code, sample_user_scores_cancer, "UserVep", None, None,
        False, qry, list_variants=True)
    assert metrics.num_variants_included == 1500


def test_compute_metrics_query_params_exclude_variant_ids(
        ve_analyzer: VEAnalyzer,
        ve_bm_query_mgr: VEBenchmarkQueryMgr,
        sample_user_scores_cancer):
    task_code = "CANCER"
    sample_scores = sample_user_scores_cancer
    qry = VEQueryCriteria(variant_ids=sample_scores[:1500],
                          include_variant_ids=False)
    metrics = ve_analyzer.compute_metrics(
        task_code, sample_scores, "UserVep", None, None,
        False, qry, list_variants=True)
    assert metrics.num_variants_included == 500


def test_compute_metrics_query_params_allele_freq(
        ve_analyzer: VEAnalyzer,
        ve_bm_query_mgr: VEBenchmarkQueryMgr,
        sample_user_scores,
        ve_reporter: VEAnalysisReporter):
    task_code = sample_user_scores[0]
    qry = VEQueryCriteria(allele_frequency_operator=">",
                          allele_frequency=1.0e-7)
    metrics = ve_analyzer.compute_metrics(
        task_code, sample_user_scores[1], "UserVep", None, None,
        False, qry, list_variants=True)
    qry1 = VEQueryCriteria(variant_ids=metrics.variants_included)
    variants = ve_bm_query_mgr.get_variants_by_task("CANCER", qry1)
    assert 0 == len(
        variants.query('ALLELE_FREQUENCY<=1.0e-7'))
    ve_reporter.write_summary(metrics)


def test_compute_metrics_query_params_multiple(
        ve_analyzer: VEAnalyzer,
        ve_bm_query_mgr: VEBenchmarkQueryMgr,
        sample_user_scores_cancer,
        ve_reporter: VEAnalysisReporter):
    task_code = "CANCER"
    sample_scores = sample_user_scores_cancer
    qry = VEQueryCriteria(allele_frequency_operator=">",
                          allele_frequency=1.0e-8,
                          variant_ids=sample_scores[:1500],
                          filter_name="Oncogene",
                          gene_symbols=['MTOR', 'PTEN'],
                          include_genes=False)
    metrics = ve_analyzer.compute_metrics(
        task_code, sample_scores, "UserVep",
        variant_effect_sources=['REVEL', 'EVE'],
        variant_query_criteria=qry, vep_min_overlap_percent=50,
        variant_vep_retention_percent=10, list_variants=True)
    ve_reporter.write_summary(metrics)


def test_compute_metrics_query_params_multiple2(
        ve_analyzer: VEAnalyzer,
        ve_bm_query_mgr: VEBenchmarkQueryMgr,
        sample_user_scores_cancer,
        ve_reporter: VEAnalysisReporter):
    task_code = "CANCER"
    sample_scores = sample_user_scores_cancer
    qry = VEQueryCriteria(
                          variant_ids=sample_scores[:1500],
                          filter_name="Oncogene",
                          gene_symbols=['MTOR', 'PTEN'],
                          include_genes=False)
    metrics = ve_analyzer.compute_metrics(
        task_code, sample_scores, "UserVep",
        variant_effect_sources=['REVEL', 'EVE'],
        variant_query_criteria=qry, vep_min_overlap_percent=50,
        variant_vep_retention_percent=10, list_variants=True)
    ve_reporter.write_summary(metrics)


def test_compute_metrics_query_params_user_vep_only(
        ve_analyzer: VEAnalyzer,
        ve_bm_query_mgr: VEBenchmarkQueryMgr,
        sample_user_scores_cancer,
        ve_reporter: VEAnalysisReporter):
    task_code = "CANCER"
    sample_scores = sample_user_scores_cancer
    qry = VEQueryCriteria(
                          variant_ids=sample_scores[:1500],
                          filter_name="Oncogene",
                          gene_symbols=['MTOR'],
                          include_genes=False)
    metrics = ve_analyzer.compute_metrics(
        task_code, sample_scores, "UserVep",
        variant_query_criteria=qry, include_variant_effect_sources=False,
        list_variants=True)
    ve_reporter.write_summary(metrics)


def test_compute_metrics_query_params_sys_veps_only(
        ve_analyzer: VEAnalyzer,
        ve_bm_query_mgr: VEBenchmarkQueryMgr,
        sample_user_scores_cancer,
        ve_reporter: VEAnalysisReporter):
    task_code = "CANCER"
    sample_scores = sample_user_scores_cancer
    qry = VEQueryCriteria(
                          variant_ids=sample_scores[:1500],
                          filter_name="Oncogene",
                          gene_symbols=['MTOR'],
                          include_genes=False)
    metrics = ve_analyzer.compute_metrics(
        task_code,
        variant_query_criteria=qry, vep_min_overlap_percent=50,
        variant_vep_retention_percent=10, list_variants=True)
    ve_reporter.write_summary(metrics)


def test_compute_metrics_sys_veps_only(
        ve_analyzer: VEAnalyzer,
        ve_bm_query_mgr: VEBenchmarkQueryMgr,
        ve_reporter: VEAnalysisReporter):
    metrics = ve_analyzer.compute_metrics(
        "CANCER",
        vep_min_overlap_percent=50,
        variant_vep_retention_percent=10, list_variants=True)
    ve_reporter.write_summary(metrics)


