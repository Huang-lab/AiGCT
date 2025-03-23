"""
Classes and methods to query the variant repository.
This layer provides an abstraction layer that sits on top of the
data access layer in the repository module.
It uses the repository module to access the raw data and includes
methods to optionally transform the data to make it more meaningful
or presentation ready to the caller.
"""

from .repository import (
    VariantEffectLabelRepository,
    VariantRepository,
    VariantTaskRepository,
    VariantEffectSourceRepository,
    VariantEffectScoreRepository,
    VariantFilterRepository,
    VARIANT_PK_COLUMNS
)
from .model import (
    VEQueryCriteria,
    VariantFilter
)
from .query_util import validate_query_criteria

import pandas as pd


def cleanup_variant_query_params(params: VEQueryCriteria):
    if params.gene_symbols is not None:
        if len(params.gene_symbols) == 0:
            params.gene_symbols = None
        else:
            if not isinstance(params.gene_symbols, pd.DataFrame):
                params.gene_symbols = pd.DataFrame(
                    {"GENE_SYMBOL": params.gene_symbols})
            if params.include_genes is None:
                params.include_genes = True
    if params.variant_ids is not None:
        if len(params.variant_ids) == 0:
            params.variant_ids = None
        else:
            if params.include_variant_ids is None:
                params.include_genes = True
    if params.column_name_map is not None:
        if len(params.column_name) == 0:
            params.column_name_map = None
    if params.allele_frequency_operator is None:
        params.allele_frequency_operator = "="
    return params


class VEBenchmarkQueryMgr:
    """
    Methods to query the variant repository
    """

    def __init__(self,
                 variant_effect_label_repo: VariantEffectLabelRepository,
                 variant_repo: VariantRepository,
                 variant_task_repo: VariantTaskRepository,
                 variant_effect_source_repo: VariantEffectSourceRepository,
                 variant_effect_score_repo: VariantEffectScoreRepository,
                 variant_filter_repo: VariantFilterRepository
                 ):
        self._variant_effect_label_repo = variant_effect_label_repo
        self._variant_repo = variant_repo
        self._variant_task_repo = variant_task_repo
        self._variant_effect_source_repo = variant_effect_source_repo
        self._variant_effect_score_repo = variant_effect_score_repo
        self._variant_filter_repo = variant_filter_repo

    def get_tasks(self) -> pd.DataFrame:
        """Get all tasks"""

        return self._variant_task_repo.get_all()

    def get_all_variants(self) -> pd.DataFrame:
        return self._variant_repo.get_all()

    def get_variants(self, qry: VEQueryCriteria) -> pd.DataFrame:
        """
        Fetch variants based on query criteria.

        Parameters
        ----------
        qry : VEQueryCriteria
            See description of VEQueryCriteria in model package.
            Specifies criteria that would limit the set of variants
            to be retrieved. The filter_name attribute is ignored.

        Returns
        -------
        DataFrame
        """
        validate_query_criteria(qry)
        return self._variant_repo.get(qry)

    def get_variant_effect_sources(self, task_code: str) -> pd.DataFrame:
        return self._variant_effect_source_repo.get_by_task(task_code)

    @staticmethod
    def _compute_variant_counts(group) -> pd.Series:
        return pd.Series(
            {"NUM_VARIANTS": len(group),
             "NUM_POSITIVE_LABELS": group["BINARY_LABEL"].sum(),
             "NUM_NEGATIVE_LABELS": (group["BINARY_LABEL"] ^ 1).sum(),
             "NUM_GENES": group["GENE_SYMBOL"].nunique()
             })

    def get_variant_effect_source_stats(
            self, task_code: str, variant_effect_sources=None,
            include_variant_effect_sources: bool = True,
            qry: VEQueryCriteria = None) -> pd.DataFrame:
        """
        Get all variant effect sources for a task along with the
        number of variants, number of positive labels,
        number of negative labels, number of genes for each source.

        Parameters
        ----------
        task_code : str

        variant_effect_sources : list, optional
            If specified it would restrict the results based on
            system supplied vep's in this list.
        include_variant_effect_sources : bool, optional
            If variant_effect_source is specified, indicates whether to
            limit the results to sources in variant_effect_sources or
            not in variant_effect_sources.
        qry : VEQueryCriteria, optional
            See description of VEQueryCriteria in model package.
            Specifies criteria that would limit the set of variants
            to be retrieved.

        Returns
        -------
        DataFrame
        """

        validate_query_criteria(qry)
        variant_labels = self._variant_effect_label_repo.get(
            task_code, qry)
        scores = self._variant_effect_score_repo.get(
            task_code, variant_effect_sources,
            include_variant_effect_sources, qry)
        scores_labels = scores.merge(variant_labels, how="inner",
                                     on=VARIANT_PK_COLUMNS)
        grouped_scores = scores_labels.groupby("SCORE_SOURCE")
        return grouped_scores.apply(
            self._compute_variant_counts,
            include_groups=False).reset_index()

    def get_all_variant_effect_source_stats(self) -> pd.DataFrame:
        scores = self._variant_effect_score_repo.get_all_for_all_tasks()
        grouped_scores = scores.groupby(["TASK_CODE"])
        return grouped_scores.agg(
            NUM_SCORE_SOURCES=('SCORE_SOURCE', 'nunique')).reset_index()

    def get_all_task_variant_effect_label_stats(self) -> pd.DataFrame:
        """
        Returns one row per task with number of variants,
        number of positive labels, number of negative labels,
        number of genes.

        Returns
        -------
        DataFrame
        """

        labels = self._variant_effect_label_repo.get_all_for_all_tasks()
        grouped_labels = labels.groupby(["TASK_CODE", "TASK_NAME"])
        return grouped_labels.apply(
            self._compute_variant_counts,
            include_groups=False).reset_index()

    def get_variant_effect_scores(self, task_code: str,
                                  variant_effect_sources=None,
                                  include_variant_effect_sources: bool = None,
                                  qry: VEQueryCriteria = None) -> pd.DataFrame:
        """
        Fetches variant effect scores for variant effect sources.

        Parameters
        ----------
        task_code : str
            task code
        variant_effect_sources : list, optional
            If specified it would restrict the results based on
            system supplied vep's in this list.
        include_variant_effect_sources : bool, optional
            If variant_effect_source is specified, indicates whether to
            limit the results to sources in variant_effect_sources or
            not in variant_effect_sources.
        qry : VEQueryCriteria, optional
            See description of VEQueryCriteria in model package.
            Specifies criteria that would limit the set of variants
            to be retrieved.

        Returns
        -------
        DataFrame
        """

        validate_query_criteria(qry)
        return self._variant_effect_score_repo.get(
            task_code,
            variant_effect_sources,
            include_variant_effect_sources,
            qry)

    def get_variants_by_task(self, task_code: str,
                             qry: VEQueryCriteria = None
                             ) -> pd.DataFrame:
        """
        Fetches variants by task. The optional parameters are 
        filter criteria used to limit the set of variants returned.

        Parameters
        ----------
        task_code : str

        qry : VEQueryCriteria, optional
            See description of VEQueryCriteria in model package.
            Specifies criteria that would limit the set of variants
            to be retrieved.

        Returns
        -------
        DataFrame
        """

        validate_query_criteria(qry)
        return self._variant_effect_label_repo.get(task_code,
                                                   qry)

    def get_variant_distribution(self, task_code: str,
                                 by: str = "gene",
                                 qry: VEQueryCriteria = None
                                 ) -> pd.DataFrame:
        """
        Fetches the distribution of variants by gene or chromsome.
        For each gene/chromosome lists number of variants for which we have
        labels along with the number of positive and negative label counts.

        Parameters
        ----------
        task_code : str
            Task code
        by : str
            Values are gene or chromosome. Specifies the type of distribution
            to return.
        qry : VEQueryCriteria, optional
            See description of VEQueryCriteria in model package.
            Specifies criteria that would limit the set of variants
            to be retrieved.

        Returns
        -------
        DataFrame
        """

        validate_query_criteria(qry)
        label_df = self._variant_effect_label_repo.get(task_code, qry)[
            ['CHROMOSOME', 'GENE_SYMBOL', 'BINARY_LABEL']]
        label_df[['POSITIVE_LABEL', 'NEGATIVE_LABEL']] = label_df.apply(
            lambda row: [row['BINARY_LABEL'], 1 ^ row['BINARY_LABEL']],
            result_type="expand", axis=1
        )
        if by == 'chromosome':
            grouped = label_df.groupby('CHROMOSOME')
        else:
            grouped = label_df.groupby('GENE_SYMBOL')
        return grouped.agg(
                NUM_POSITIVE_LABELS=pd.NamedAgg(column='POSITIVE_LABEL',
                                                aggfunc='sum'),
                NUM_NEGATIVE_LABELS=pd.NamedAgg(column='NEGATIVE_LABEL',
                                                aggfunc='sum')
                            ).reset_index()

    def get_variant_filter(
            self, task_code: str, filter_name: str) -> VariantFilter:
        """
        Return a variant filter for a task by name.

        Returns
        -------
        VariantFilter
            Object containing list of genes/variant id's included in the
            filter. See description of the object.
        """
        return self._variant_filter_repo.get_by_task_filter_name(
            task_code, filter_name
        )

    def get_all_variant_filters(self, task_code: str) ->\
            dict[str, pd.DataFrame]:
        """
        Return basic descriptive information about all variant filters
        for a task.

        Returns
        -------
        dict[str, pd.DataFrame]
            A dictionary of 3 data frames with the following keys:
            filter_df - Data frame of filters containing CODE, NAME,
            DESCRIPTION, etc.
            filter_gene_df - Data frame of genes associated with each
            filter
            filter_variant_df - Data frame of variants associated with
            each filter
        """
        return self._variant_filter_repo.get_by_task(task_code)
