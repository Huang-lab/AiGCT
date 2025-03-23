import context  # noqa: F401
from aigct.model import VEQueryCriteria
from aigct.query import VEBenchmarkQueryMgr


def test_query_criteria(ve_bm_query_mgr: VEBenchmarkQueryMgr,
                        sample_user_scores_cancer):
    qry = VEQueryCriteria(
                          variant_ids=sample_user_scores_cancer[:1500],
                          filter_name="Oncogene",
                          gene_symbols=['MTOR', 'PTEN'],
                          include_genes=False)
    variants_df = ve_bm_query_mgr.get_all_variants()
    assert len(variants_df) > 0
    tasks_df = ve_bm_query_mgr.get_tasks()
    assert len(tasks_df) > 0
    sources_df = ve_bm_query_mgr.get_variant_effect_sources("CANCER")
    assert len(sources_df) > 0
    score_df = ve_bm_query_mgr.get_variant_effect_scores("CANCER", qry=qry)
    assert len(score_df) > 0
    filter_dfs = ve_bm_query_mgr.get_variant_filter("CANCER", "Oncogene")
    assert len(filter_dfs.filter) > 0 and (
        len(filter_dfs.filter_genes) > 0 or
        len(filter_dfs.filter_variants) > 0)
    filters = ve_bm_query_mgr.get_all_variant_filters("CANCER")
    assert (len(filters["filter_df"]) > 0 and
            (len(filters["filter_gene_df"]) > 0 or
            len(filters["filter_variant_df"]) > 0))
    variant_df = ve_bm_query_mgr.get_variants(qry)
    assert len(variant_df) > 0
    variants_by_task_df = ve_bm_query_mgr.get_variants_by_task("CANCER", qry)
    assert (len(variants_by_task_df) > 0 and len(variants_by_task_df) <
            len(variants_df))


def test_query_source_stats(ve_bm_query_mgr: VEBenchmarkQueryMgr,
                            sample_user_scores_cancer):
    qry = VEQueryCriteria(
                          variant_ids=sample_user_scores_cancer[:1500],
                          filter_name="Oncogene",
                          gene_symbols=['MTOR', 'PTEN'],
                          include_genes=False)
    stats = ve_bm_query_mgr.get_variant_effect_source_stats(
        "CANCER", qry=qry)
    assert len(stats) > 0
    stats = ve_bm_query_mgr.get_variant_effect_source_stats(
        "CANCER", ["ALPHAM", "REVEL", "ESMB1"])
    assert len(stats) > 0


def test_query_dist(ve_bm_query_mgr: VEBenchmarkQueryMgr,
                    sample_user_scores_cancer):
    qry = VEQueryCriteria(
                          variant_ids=sample_user_scores_cancer[:1500],
                          filter_name="Oncogene",
                          gene_symbols=['MTOR', 'PTEN'],
                          include_genes=False)
    dist_gene = ve_bm_query_mgr.get_variant_distribution("CANCER")
    assert len(dist_gene) > 0
    dist_chr = ve_bm_query_mgr.get_variant_distribution("CANCER",
                                                        "chromosome")
    assert len(dist_chr) > 0
    dist_qry = ve_bm_query_mgr.get_variant_distribution("CANCER",
                                                        qry=qry)
    assert len(dist_qry) > 0 and len(dist_qry) < len(dist_gene)


def test_query_get_all_labels_for_all_tasks(ve_bm_query_mgr: VEBenchmarkQueryMgr):
    labels = ve_bm_query_mgr.get_all_task_variant_effect_label_stats()
    assert len(labels) > 0


def test_get_all_variant_effect_source_stats(
        ve_bm_query_mgr: VEBenchmarkQueryMgr):
    score_stats = ve_bm_query_mgr.get_all_variant_effect_source_stats()
    assert len(score_stats) > 0