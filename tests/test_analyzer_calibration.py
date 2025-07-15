import context  # noqa: F401
import pytest
from aigct.analyzer import VEAnalyzer
from aigct.model import VEQueryCriteria
from aigct.query import VEBenchmarkQueryMgr
from aigct.reporter import VEAnalysisReporter
from aigct.pd_util import filter_dataframe_by_list
from aigct.repository import VARIANT_PK_COLUMNS

import pandas as pd

# from aigct.model import VariantId  # noqa: F401

test_tasks = ["CANCER", "CLINVAR", "DDD", "CHD", "ASD", "ADRD" ]
# test_tasks = ["CLINVAR" ]


def test_compute_calibration_metrics_basic(
        ve_analyzer: VEAnalyzer,
        sample_user_scores_cancer):
    user_scores = sample_user_scores_cancer
    metrics = ve_analyzer.compute_calibration_metrics(
        "CANCER", user_scores, "UserVep", column_name_map=None,
        variant_effect_source=None, variant_query_criteria=None,
        pathogenic_fraction_bins=10)
    metrics
