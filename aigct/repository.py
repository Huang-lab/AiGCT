"""
Data access layer methods for accessing variant repository.
Classes here provide an encapsulation layer to hide the internal
details of the repository structure.
"""

import os
import pandas as pd
import dask.dataframe as dd
from .util import ParameterizedSingleton
import threading
from dataclasses import dataclass, field
from .pd_util import (
    filter_dataframe_by_list,
    build_dataframe_where_clause,
    merge_by_chunks
)
from .model import VariantFilter, VEQueryCriteria

MERGE_CHUNK_SIZE = 500000

TASK_SUBFOLDER = {
    "CANCER": "CANCER"
}

DATA_FOLDER = "data"
TASK_FOLDERS = [os.path.join(DATA_FOLDER, task) for task in ["CANCER", "ADRD",
                                                             "CHD", "DDD"]]


@dataclass
class TableDef:
    folder: str
    file_name: str
    pk_columns: list[str]
    non_pk_columns: list[str]
    columns: list[str] = field(init=False),
    full_file_name: str = field(init=False)

    def __post_init__(self):
        self.columns = self.pk_columns + self.non_pk_columns
        self.full_file_name = os.path.join(self.folder, self.file_name)
    


VARIANT_PK_COLUMNS = [
    "GENOME_ASSEMBLY",
    "CHROMOSOME",
    "POSITION",
    "REFERENCE_NUCLEOTIDE",
    "ALTERNATE_NUCLEOTIDE"
]
VARIANT_NON_PK_COLUMNS = [
    "PRIOR_GENOME_ASSEMBLY",
    "PRIOR_CHROMOSOME",
    "PRIOR_POSITION",
    "PRIOR_PRIOR_GENOME_ASSEMBLY",
    "PRIOR_PRIOR_CHROMOSOME",
    "PRIOR_PRIOR_POSITION",
    "REFERENCE_AMINO_ACID",
    "ALTERNATE_AMINO_ACID",
    "AMINO_ACID_POSITION",
    "RS_DBSNP",
    "GENE_SYMBOL",
    "ENSEMBL_GENE_ID",
    "ENSEMBL_TRANSCRIPT_ID",
    "ENSEMBL_PROTEIN_ID",
    "ALLELE_FREQUENCY_SOURCE",
    "ALLELE_FREQUENCY"
]
VARIANT_TABLE_DEF = TableDef(DATA_FOLDER, "variant.csv",
                             VARIANT_PK_COLUMNS,
                             VARIANT_NON_PK_COLUMNS)

VARIANT_LABEL_NON_PK_COLUMNS = [
    "LABEL_SOURCE",
    "RAW_LABEL",
    "BINARY_LABEL"
]
VARIANT_EFFECT_LABEL_TABLE_DEF = TableDef(DATA_FOLDER,
                                          "variant_effect_label.csv",
                                          VARIANT_PK_COLUMNS,
                                          VARIANT_LABEL_NON_PK_COLUMNS)

VARIANT_EFFECT_SCORE_PK_COLUMNS = VARIANT_PK_COLUMNS + ["SCORE_SOURCE"]
VARIANT_EFFECT_SCORE_NON_PK_COLUMNS = [
    "RAW_SCORE",
    "RANK_SCORE"
]
VARIANT_EFFECT_SCORE_TABLE_DEF = TableDef(DATA_FOLDER,
                                          "variant_effect_score.csv",
                                          VARIANT_EFFECT_SCORE_PK_COLUMNS,
                                          VARIANT_EFFECT_SCORE_NON_PK_COLUMNS)

VARIANT_TASK_TABLE_DEF = TableDef(DATA_FOLDER,
                                  "variant_task.csv",
                                  ["CODE"], ["NAME", "DESCRIPTION"])

VARIANT_EFFECT_SOURCE_TABLE_DEF =\
    TableDef(DATA_FOLDER,
             "variant_effect_source.csv", ["CODE"],
             ["NAME", "SOURCE_TYPE", "DESCRIPTION"])

VARIANT_DATA_SOURCE_TABLE_DEF =\
    TableDef(DATA_FOLDER,
             "variant_data_source.csv", ["CODE"],
             ["NAME", "DESCRIPTION"])

VARIANT_FILTER_TABLE_DEF =\
    TableDef(DATA_FOLDER,
             "variant_filter.csv", ["CODE"],
             ["NAME", "DESCRIPTION", "INCLUDE_GENES",
              "INCLUDE_VARIANTS"])

VARIANT_FILTER_GENE_TABLE_DEF =\
    TableDef(DATA_FOLDER,
             "variant_filter_gene.csv", ["FILTER_CODE", "GENE_SYMBOL"],
             [])

VARIANT_FILTER_VARIANT_TABLE_DEF =\
    TableDef(DATA_FOLDER,
             "variant_filter_variant.csv",
             ["FILTER_CODE"] + VARIANT_PK_COLUMNS, [])

TABLE_DEFS = {
    "VARIANT_TASK": VARIANT_TASK_TABLE_DEF,
    "VARIANT_EFFECT_SOURCE": VARIANT_EFFECT_SOURCE_TABLE_DEF,
    "VARIANT": VARIANT_TABLE_DEF,
    "VARIANT_EFFECT_LABEL": VARIANT_EFFECT_LABEL_TABLE_DEF,
    "VARIANT_DATA_SOURCE": VARIANT_DATA_SOURCE_TABLE_DEF,
    "VARIANT_EFFECT_SCORE": VARIANT_EFFECT_SCORE_TABLE_DEF,
    "VARIANT_FILTER": VARIANT_FILTER_TABLE_DEF,
    "VARIANT_FILTER_GENE": VARIANT_FILTER_GENE_TABLE_DEF,
    "VARIANT_FILTER_VARIANT": VARIANT_FILTER_VARIANT_TABLE_DEF,
}


def read_repo_csv(file: str) -> pd.DataFrame:
    return pd.read_csv(file,
                       dtype={"CHROMOSOME": str, "BINARY_LABEL": int})


class RepoSessionContext:

    def __init__(self, data_folder_root: str,
                 table_defs: dict[str, TableDef]):
        self._data_folder_root = data_folder_root
        self._table_defs = table_defs

    @property
    def data_folder_root(self):
        return self._data_folder_root

    def table_def(self, table_name: str):
        return self._table_defs[table_name]

    def table_file(self, table_name: str, task: str = None):
        if task:
            return os.path.join(self._data_folder_root,
                                DATA_FOLDER, task,
                                self._table_defs[table_name].file_name)
        else:
            return os.path.join(self._data_folder_root,
                                DATA_FOLDER,
                                self._table_defs[table_name].file_name)


class VariantEffectLabelCache(ParameterizedSingleton):
    """
    Caches the variant csv file in a dataframe. Implements the singleton
    pattern to ensure there is only one instance of the cached dataframe.
    We use an _init_once method rather than the normal __init__ method
    as required by the ParameterizedSingleton class.
    """

    def _init_once(self, data_folder_root: str):
        self._data_folder_root = data_folder_root
        self._lock = threading.Lock()
        self._cache = dict()

    def get_data_frame(self, task_code: str):
        if task_code not in self._cache:
            with self._lock:
                if task_code not in self._cache:
                    self._cache[task_code] = read_repo_csv(
                        os.path.join(self._data_folder_root, DATA_FOLDER,
                                     task_code,
                                     VARIANT_EFFECT_LABEL_TABLE_DEF.file_name))
        return self._cache[task_code]


class DataCache(ParameterizedSingleton):
    """
    Caches a repository csv file in a dataframe. Implements the singleton
    pattern to ensure there is only one instance of the cached dataframe.
    We use an _init_once method rather than the normal __init__ method
    as required by the ParameterizedSingleton class.
    """

    def _init_once(self, data_folder_root: str, table_def: TableDef):
        self._data_folder_root = data_folder_root
        self._table_def = table_def
        self._lock = threading.Lock()
        self._data_frame = None

    @property
    def data_frame(self):
        if self._data_frame is None:
            with self._lock:
                if self._data_frame is None:
                    self._data_frame = read_repo_csv(
                        os.path.join(self._data_folder_root,
                                     self._table_def.full_file_name))
        return self._data_frame


class TaskBasedDataCache(ParameterizedSingleton):
    """
    Caches a repository csv file in a dataframe. Maintains a separate
    cache for each task in a dict. Implements the singleton
    pattern to ensure there is only one instance of the cached dataframe.
    We use an _init_once method rather than the normal __init__ method
    as required by the ParameterizedSingleton class.
    """

    def _init_once(self, data_folder_root: str, table_def: TableDef,
                   disable_cache: bool = False):
        self._data_folder_root = data_folder_root
        self._table_def = table_def
        self._cache = dict()
        self._lock = threading.Lock()
        self._disable_cache = disable_cache

    def get_data_frame(self, task_code: str):
        if task_code not in self._cache:
            with self._lock:
                if task_code not in self._cache:
                    data_frame = read_repo_csv(
                        os.path.join(self._data_folder_root, DATA_FOLDER,
                                     task_code,
                                     self._table_def.file_name))
                    if self._disable_cache:
                        return data_frame
                    self._cache[task_code] = data_frame
        return self._cache[task_code]


class TaskDataCache(DataCache):
    """
    Caches the variant csv file in a dataframe. Implements the singleton
    pattern to ensure there is only one instance of the cached dataframe.
    We use an _init_once method rather than the normal __init__ method
    as required by the ParameterizedSingleton class.
    """

    def _init_once(self, data_folder_root: str):
        super()._init_once(data_folder_root, VARIANT_TASK_TABLE_DEF)


class VariantEffectScoreCache(TaskBasedDataCache):
    """
    Caches the variant csv file in a dataframe. Implements the singleton
    pattern to ensure there is only one instance of the cached dataframe.
    We use an _init_once method rather than the normal __init__ method
    as required by the ParameterizedSingleton class.
    """

    def _init_once(self, data_folder_root: str, disable_cache: bool = False):
        super()._init_once(data_folder_root, VARIANT_EFFECT_SCORE_TABLE_DEF,
                           disable_cache)


class VariantCache(DataCache):
    """
    Caches the variant csv file in a dataframe. Implements the singleton
    pattern to ensure there is only one instance of the cached dataframe.
    We use an _init_once method rather than the normal __init__ method
    as required by the ParameterizedSingleton class.
    """

    def _init_once(self, data_folder_root: str):
        super()._init_once(data_folder_root, VARIANT_TABLE_DEF)


class VariantTaskCache(DataCache):
    """
    Caches the variant csv file in a dataframe. Implements the singleton
    pattern to ensure there is only one instance of the cached dataframe.
    We use an _init_once method rather than the normal __init__ method
    as required by the ParameterizedSingleton class.
    """

    def _init_once(self, data_folder_root: str):
        super()._init_once(data_folder_root, VARIANT_TASK_TABLE_DEF)


class VariantEffectSourceCache(DataCache):

    def _init_once(self, data_folder_root: str):
        super()._init_once(data_folder_root, VARIANT_EFFECT_SOURCE_TABLE_DEF)


class VariantFilterCache(ParameterizedSingleton):

    def _init_once(self, data_folder_root: str):
        self._data_folder_root = data_folder_root
        self._lock = threading.Lock()
        self._cache = dict()

    def get_data_frames(self, task_code: str) -> dict:
        if task_code not in self._cache:
            with self._lock:
                if task_code not in self._cache:
                    folder = os.path.join(self._data_folder_root, DATA_FOLDER,
                                          task_code)
                    cache_dict = dict()
                    cache_dict["filter_df"] = read_repo_csv(
                        os.path.join(folder,
                                     VARIANT_FILTER_TABLE_DEF.file_name))
                    cache_dict["filter_gene_df"] = read_repo_csv(
                        os.path.join(folder,
                                     VARIANT_FILTER_GENE_TABLE_DEF.file_name))
                    cache_dict["filter_variant_df"] = read_repo_csv(
                        os.path.join(
                            folder,
                            VARIANT_FILTER_VARIANT_TABLE_DEF.file_name))

                    self._cache[task_code] = cache_dict
        return self._cache[task_code]


class VariantEffectSourceRepository:

    def __init__(self, session_context: RepoSessionContext,
                 variant_effect_score_repo):
        self._cache = VariantEffectSourceCache(
            session_context.data_folder_root)
        self._variant_effect_score_repo = variant_effect_score_repo

    def get_all(self) -> pd.DataFrame:
        return self._cache.data_frame.copy(deep=True)

    def get_by_task(self, task_code: str) -> pd.DataFrame:
        score_sources = self._variant_effect_score_repo.get_all_by_task(
            task_code)['SCORE_SOURCE'].unique()
        source_df = self._cache.data_frame.copy(deep=True)
        return source_df[source_df['CODE'].isin(score_sources)]

    def get_by_code(self, codes: list[str]) -> pd.DataFrame:
        return filter_dataframe_by_list(
            self._cache.data_frame, codes,
            VARIANT_EFFECT_SOURCE_TABLE_DEF.pk_columns)


class VariantTaskRepository:

    def __init__(self, session_context: RepoSessionContext):
        self._cache = VariantTaskCache(session_context.data_folder_root)

    def get_all(self) -> pd.DataFrame:
        return self._cache.data_frame.copy(deep=True)


class VariantFilterRepository:

    def __init__(self, session_context: RepoSessionContext):
        self._cache = VariantFilterCache(session_context.data_folder_root)

    def get_by_task(self, task_code: str) -> dict[str, pd.DataFrame]:
        return self._cache.get_data_frames(task_code)

    def get_by_task_filter_name(
            self, task_code: str, filter_name: str) -> VariantFilter:
        filter_dfs = self._cache.get_data_frames(task_code)
        filter = filter_dfs["filter_df"].query(f"NAME == '{filter_name}'")
        if len(filter) == 0:
            return None
        filter = filter.iloc[0]
        filter_code = filter["CODE"]
        filter_genes = filter_dfs["filter_gene_df"].query(
            f"FILTER_CODE == '{filter_code}'")
        filter_variants = filter_dfs["filter_variant_df"].query(
            f"FILTER_CODE == '{filter_code}'")
        if len(filter_genes) == 0 and len(filter_variants) == 0:
            return None
        return VariantFilter(filter, filter_genes, filter_variants)

    def get_by_task_filter_names(
            self, task_code: str, filter_names: list[str]) -> list[
                VariantFilter]:
        variant_filters = []
        for filter_name in filter_names:
            variant_filter = self.get_by_task_filter_name(task_code, filter_name)
            if variant_filter is None:
                raise Exception(f"Invalid filter name: {filter_name}")
            variant_filters.append(variant_filter)
        return variant_filters


def query_by_filter(query_df: pd.DataFrame,
                    filter: pd.Series,
                    filter_gene_df: pd.DataFrame,
                    filter_variant_df: pd.DataFrame) -> pd.DataFrame:

    if len(filter_gene_df) > 0:
        query_df = filter_dataframe_by_list(query_df, filter_gene_df,
                                            ['GENE_SYMBOL'], None,
                                            filter["INCLUDE_GENES"] ==
                                            'Y')
    if len(filter_variant_df) > 0:
        query_df = filter_dataframe_by_list(query_df, filter_variant_df,
                                            VARIANT_PK_COLUMNS, None,
                                            filter["INCLUDE_VARIANTS"] == 'Y')
    return query_df


def query_by_filters(query_df: pd.DataFrame, filters: list[VariantFilter]
                     ) -> pd.DataFrame:
    results = []
    for filter in filters:
        result = query_by_filter(query_df, filter.filter,
                                 filter.filter_genes,
                                 filter.filter_variants)
        results.append(result)
    return pd.concat(results).drop_duplicates()


class VariantRepository:

    def __init__(self, session_context: RepoSessionContext):
        self._cache = VariantCache(session_context.data_folder_root)

    def get_all(self) -> pd.DataFrame:
        return self._cache.data_frame.copy(deep=True)

    def get(self, qry: VEQueryCriteria) -> pd.DataFrame:
        """
        Fetches variants. The optional parameters are filter criteria used to
        limit the set of variants returned.
        """
        where_clause = build_dataframe_where_clause(
            {"ALLELE_FREQUENCY": [qry.allele_frequency_operator,
                                  qry.allele_frequency]})
        variant_df = self._cache.data_frame
        if where_clause != "":
            variant_df = variant_df.query(where_clause)
        if qry.gene_symbols is not None:
            variant_df = filter_dataframe_by_list(variant_df, qry.gene_symbols,
                                                  ['GENE_SYMBOL'],
                                                  in_list=qry.include_genes)
        if qry.variant_ids is not None:
            variant_df = filter_dataframe_by_list(variant_df, qry.variant_ids,
                                                  VARIANT_PK_COLUMNS,
                                                  qry.column_name_map,
                                                  qry.include_variant_ids)
        return variant_df


class VariantEffectLabelRepository:

    def __init__(self, session_context: RepoSessionContext,
                 variant_task_repo: VariantTaskRepository,
                 variant_repo: VariantRepository,
                 filter_repo: VariantFilterRepository):
        self._cache = VariantEffectLabelCache(session_context.data_folder_root)
        self._task_repo = variant_task_repo
        self._filter_repo = filter_repo
        self._variant_repo = variant_repo

    def get_all_by_task(self, task_code: str) -> pd.DataFrame:
        label_df = self._cache.get_data_frame(task_code)
        variant_df = self._variant_repo.get(
            VEQueryCriteria(variant_ids=label_df[VARIANT_PK_COLUMNS]))
        return label_df.merge(variant_df, on=VARIANT_PK_COLUMNS,
                              how="inner")

    def get_all_for_all_tasks(self) -> pd.DataFrame:
        tasks_df = self._task_repo.get_all()
        labels_dfs = []
        for row in tasks_df.itertuples():
            labels_df = self.get_all_by_task(row.CODE)
            labels_df["TASK_CODE"] = row.CODE
            labels_df["TASK_NAME"] = row.NAME
            labels_dfs.append(labels_df)
        return pd.concat(labels_dfs)

    def get(self, task_code: str,
            qry: VEQueryCriteria = None) -> pd.DataFrame:
        """
        Fetches variant effect labels.
        """
        label_df = self._cache.get_data_frame(task_code)
        if qry is not None:
            if qry.variant_ids is not None and len(qry.variant_ids) > 0:
                label_df = filter_dataframe_by_list(
                    label_df, qry.variant_ids,
                    VARIANT_PK_COLUMNS,
                    in_list=qry.include_variant_ids)
            variant_df = self._variant_repo.get(qry)
        else:
            variant_df = self._variant_repo.get_all()
        merge_df = label_df.merge(variant_df, how="inner",
                                  on=VARIANT_PK_COLUMNS)
        if qry is not None and qry.filter_names is not None:
            if isinstance(qry.filter_names, str):
                filter_names = [qry.filter_names]
            else:
                filter_names = qry.filter_names
            filters = self._filter_repo.get_by_task_filter_names(
                task_code, filter_names)
            merge_df = query_by_filters(merge_df, filters)
        return merge_df


class VariantEffectScoreRepository:

    def __init__(self, session_context: RepoSessionContext,
                 task_repo: VariantTaskRepository,
                 variant_repo: VariantRepository,
                 filter_repo: VariantFilterRepository):
        self._cache = VariantEffectScoreCache(session_context.data_folder_root,
                                              disable_cache=True)
        self._task_repo = task_repo
        self._filter_repo = filter_repo
        self._variant_repo = variant_repo

    def get_all_by_task(self, task_code: str) -> pd.DataFrame:
        score_df = self._cache.get_data_frame(task_code)
        query_criteria = VEQueryCriteria(variant_ids=score_df[
            VARIANT_PK_COLUMNS].drop_duplicates(),
            include_variant_ids=True)
        variant_df = self._variant_repo.get(query_criteria)
        return merge_by_chunks(score_df, variant_df,
                               chunk_size=MERGE_CHUNK_SIZE,
                               on=VARIANT_PK_COLUMNS,
                               how="inner")

    def get_all_by_task_slim(self, task_code: str) -> pd.DataFrame:
        score_df = self._cache.get_data_frame(task_code)
        return score_df

    def get(self, task_code: str,
            variant_effect_sources: list[str] | str = None,
            include_variant_effect_sources: bool = True,
            qry: VEQueryCriteria = None) -> pd.DataFrame:

        score_df = self._cache.get_data_frame(task_code)
        if (variant_effect_sources is not None and
                len(variant_effect_sources) > 0):
            score_df = filter_dataframe_by_list(score_df,
                                                variant_effect_sources,
                                                "SCORE_SOURCE",
                                                None,
                                                include_variant_effect_sources)
        if qry is None:
            return score_df
        if qry.variant_ids is not None:
            score_df = filter_dataframe_by_list(score_df, qry.variant_ids,
                                                VARIANT_PK_COLUMNS,
                                                qry.column_name_map,
                                                qry.include_variant_ids)
        variant_df = self._variant_repo.get(qry)
        merge_df = merge_by_chunks(score_df, variant_df,
                                   chunk_size=MERGE_CHUNK_SIZE,
                                   how="inner",
                                   on=VARIANT_PK_COLUMNS)
        if qry.filter_names is not None:
            if isinstance(qry.filter_names, str):
                filter_names = [qry.filter_names]
            else:
                filter_names = qry.filter_names
            filters = self._filter_repo.get_by_task_filter_names(
                task_code, filter_names)
            merge_df = query_by_filters(merge_df, filters)
        return merge_df[VARIANT_EFFECT_SCORE_TABLE_DEF.columns]

    def get_all_for_all_tasks(self) -> pd.DataFrame:
        tasks_df = self._task_repo.get_all()
        scores_dfs = []
        for row in tasks_df.itertuples():
            scores_df = self.get_all_by_task_slim(row.CODE)
            scores_df["TASK_CODE"] = row.CODE
            scores_df["TASK_NAME"] = row.NAME
            scores_dfs.append(scores_df)
        return pd.concat(scores_dfs)


