from dataclasses import dataclass
import pandas as pd
from .model import (
    PkViolations,
    TaskPkViolations
)

from .repository import (
    VariantEffectLabelRepository,
    VariantRepository,
    VariantTaskRepository,
    VariantEffectSourceRepository,
    VariantEffectScoreRepository,
    VariantFilterRepository,
    VARIANT_PK_COLUMNS,
    VARIANT_EFFECT_SOURCE_TABLE_DEF,
    VARIANT_EFFECT_LABEL_TABLE_DEF,
    VARIANT_EFFECT_SCORE_PK_COLUMNS,
    VARIANT_FILTER_GENE_TABLE_DEF,
    VARIANT_FILTER_VARIANT_TABLE_DEF,
    VARIANT_FILTER_TABLE_DEF,
)


class VEDataValidator:

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

    def _validate_pk(self, df: pd.DataFrame, pk_columns: list[str]
                     ) -> pd.DataFrame:
        grouped_df = df.groupby(pk_columns).size()
        dups = grouped_df[grouped_df > 1]
        return dups.reset_index()

    def validate_pks(self):
        variants = self._variant_repo.get_all()
        variant_dups = self._validate_pk(variants, VARIANT_PK_COLUMNS)
        ve_sources = self._variant_effect_source_repo.get_all()
        ve_source_dups = self._validate_pk(
            ve_sources, VARIANT_EFFECT_SOURCE_TABLE_DEF.pk_columns)
        tasks = self._variant_task_repo.get_all()
        task_violations: dict[str, TaskPkViolations] = {}
        task_dups_found = False
        for i, task in tasks.iterrows():
            task_code = task['CODE']
            labels = self._variant_effect_label_repo.get_all_by_task(
                task_code
            )
            label_dups = self._validate_pk(
                labels, VARIANT_EFFECT_LABEL_TABLE_DEF.pk_columns)
            scores = self._variant_effect_score_repo.get_all_by_task(
                    task_code
            )
            score_dups = self._validate_pk(
                scores, VARIANT_EFFECT_SCORE_PK_COLUMNS)
            filters = self._variant_filter_repo.get_by_task(task_code)
            filter_dups = self._validate_pk(
                    filters["filter_df"],
                    VARIANT_FILTER_TABLE_DEF.pk_columns)
            filter_gene_dups = self._validate_pk(
                    filters["filter_gene_df"],
                    VARIANT_FILTER_GENE_TABLE_DEF.pk_columns)
            filter_variant_dups = self._validate_pk(
                    filters["filter_variant_df"],
                    VARIANT_FILTER_VARIANT_TABLE_DEF.pk_columns)
            dups_found = bool(len(label_dups) | len(score_dups) |
                              len(filter_dups) | len(filter_gene_dups) |
                              len(filter_variant_dups))
            if dups_found:
                task_dups_found = True
                task_violations[task_code] = TaskPkViolations(
                    dups_found,
                    label_dups, score_dups, filter_dups,
                    filter_gene_dups, filter_variant_dups)
        dups_found = bool(task_dups_found | len(variant_dups) |
                          len(ve_source_dups))
        return PkViolations(dups_found, variant_dups, ve_source_dups,
                            task_violations)
