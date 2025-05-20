import pytest
import random
import context  # noqa: F401
from aigct.container import VEBenchmarkContainer
from aigct.repository import (
    VARIANT_PK_COLUMNS
)
TEST_TASK_VEP_CODES = [("CANCER","REVEL"), 
                  ("CLINVAR","VAR_R"),
                  ("DDD","VAR_R"),
                  ("CHD","VAR_R"),
                  ("ASD","VAR_R"),
                   ("ADRD","VAR_R")]


def get_ve_bm_container():
    return VEBenchmarkContainer()


def generate_sample_user_scores(task_code, vep_code):
    user_scores_df = get_ve_bm_container()._score_repo.get(
        task_code, vep_code)
    random_idxs = random.sample(list(range(len(user_scores_df))),
                                min(len(user_scores_df), 2000))
    user_variants = user_scores_df.iloc[random_idxs].copy()
    user_variants['RANK_SCORE'] = user_variants['RANK_SCORE'].apply(
        lambda scor: scor + (random.uniform(0.05, 0.15) *
                             random.sample([1, -1], 1)[0])
        if scor < 0.84 and scor > 0.16 else scor)
    # user_variants["CHROMOSOME"] = user_variants["CHROMOSOME"].astype(str)
    return user_variants


@pytest.fixture
def ve_bm_container():
    return get_ve_bm_container()


@pytest.fixture
def ve_analyzer(ve_bm_container):
    return ve_bm_container.analyzer


@pytest.fixture
def ve_reporter(ve_bm_container):
    return ve_bm_container.reporter


@pytest.fixture
def ve_plotter(ve_bm_container):
    return ve_bm_container.plotter


def generate_random_floats(n, start, end):
    return [random.uniform(start, end) for _ in range(n)]

@pytest.fixture(scope="module", params=TEST_TASK_VEP_CODES)
def task_vep_code(request):
    return request.param


@pytest.fixture
def sample_user_scores(ve_bm_container, task_vep_code):
    task_code = task_vep_code[0]
    vep_code = task_vep_code[1]
    user_variants = generate_sample_user_scores(task_code, vep_code)
    return task_code, user_variants


@pytest.fixture
def sample_user_scores_cancer(ve_bm_container):
    return generate_sample_user_scores("CANCER", "REVEL")


@pytest.fixture
def sample_user_scores_col_name_map(ve_bm_container, sample_user_scores):
    user_scores_df = ve_bm_container._score_repo.get("CANCER",
                                                          "REVEL")
    user_scores_df["RANK_SCORE"] = user_scores_df["RANK_SCORE"].apply(
        lambda sc: sc + 0.1 if sc < 0.9 else 0.95)
    col_name_map = {"CHROMOSOME": "chr", "POSITION": "pos",
                    "REFERENCE_NUCLEOTIDE": "ref",
                    "ALTERNATE_NUCLEOTIDE": "alt",
                    "GENOME_ASSEMBLY": "assembly",
                    "RANK_SCORE": "score"}
    user_scores_df.rename(columns=col_name_map, inplace=True)
    col_name_map = {val: key for key, val in col_name_map.items()}
    return user_scores_df, col_name_map


@pytest.fixture
def ve_bm_query_mgr(ve_bm_container):
    return ve_bm_container.query_mgr


@pytest.fixture
def ve_analysis_exporter(ve_bm_container):
    return ve_bm_container.exporter


@pytest.fixture
def ve_data_validator(ve_bm_container):
    return ve_bm_container.data_validator
 


@pytest.fixture
def sample_veps(ve_bm_container):
    srcs = ve_bm_container.query_mgr.get_variant_effect_sources()
    return list(srcs[:5]["CODE"])

