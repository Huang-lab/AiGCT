Database Query Guide
====================

Do the following to prepare for querying the database. <config> is the directory
where the config file is stored as specified during the installation::

    import pandas as pd
    from aigct.container import VEBenchmarkContainer

    container = VEBenchmarkContainer("<config>/aigct.yaml")

    query_mgr = container.query_mgr

query_mgr is an instance of the aigct.query.VEBenchmarkQueryMgr class that contains
the query methods. Most methods return the results as a dataframe. Here are some
some of the methods available::
    
    # Fetch all tasks (categories of variants that we have data for)

    tasks_df = query_mgr.get_tasks()

    # Fetch all variants across all tasks

    variants_df = query_mgr.get_all_variants()

    # Fetch all variant effect sources ,i.e. VEP's for which we have scores,
    # for the CANCER task

    veps_df = query_mgr.get_variant_effect_sources("CANCER")

    # Fetch statics for a specific set of variant effect sources
    # for the CANCER task

    vep_stats_df = query_mgr.get_variant_effect_source_stats("CANCER",
        ["ALPHAM", "REVEL", "EVE"])

Below we illustrate how detailed selection criteria can be specified for a query method::

    from aigct.model import VEQueryCriteria

    # Fetch scores for all variant effect sources for the CANCER task.
    # Limit to variants found in the MTOR, PLCH2, PIK3CD genes.
    # VEQueryCriteria allows you to specify detailed selection criteria
    # for the query.

    selection_criteria = VeQueryCriteria(gene_symbols=["MTOR", "PLCH2", "PIK3CD"])

    scores_df = query_mgr.get_variant_effect_scores("CANCER",
        qry=selection_criteria)

Detailed information about all of the query methods available can be found in the
API Documentation for the aigct.query.VEBenchmarkQueryMgr class.

