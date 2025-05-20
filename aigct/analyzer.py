import pandas as pd
import numpy as np
from scipy import stats
from pandas.api.types import is_string_dtype
from sklearn.metrics import roc_curve, roc_auc_score, precision_recall_curve, auc
from .model import VEQueryCriteria, VEAnalysisResult
from .repository import (
    VariantEffectScoreRepository,
    VariantEffectLabelRepository,
    VariantEffectSourceRepository,
    VARIANT_PK_COLUMNS,
)
from .pd_util import filter_dataframe_by_list
from .query_util import validate_query_criteria

VARIANT_EFFECT_SCORE_COLS = ["SCORE_SOURCE"] + VARIANT_PK_COLUMNS + ["RANK_SCORE"]


class VEAnalyzer:

    def __init__(
        self,
        variant_effect_score_repo: VariantEffectScoreRepository,
        variant_effect_label_repo: VariantEffectLabelRepository,
        variant_effect_source_repo: VariantEffectSourceRepository,
    ):
        self._variant_effect_score_repo = variant_effect_score_repo
        self._variant_effect_label_repo = variant_effect_label_repo
        self._variant_effect_source_repo = variant_effect_source_repo

    def validate_user_ve_scores(self, user_ve_scores: pd.DataFrame):
        if (
            user_ve_scores is not None
            and len(user_ve_scores) > 0
            and not is_string_dtype(user_ve_scores["CHROMOSOME"])
        ):
            raise Exception("CHROMOSOME column values must all be string type")

    def get_analysis_scores_and_labels(
        self,
        task_code: str,
        user_ve_scores: pd.DataFrame = None,
        user_vep_name: str = "USER",
        column_name_map: dict = None,
        variant_effect_sources: list[str] = None,
        include_variant_effect_sources: bool = None,
        variant_query_criteria: VEQueryCriteria = None,
        vep_min_overlap_percent: float = 0,
        variant_vep_retention_percent: float = 0,
    ) -> pd.DataFrame:

        if vep_min_overlap_percent is None:
            vep_min_overlap_percent = 0
        if variant_vep_retention_percent is None:
            variant_vep_retention_percent = 0
        if user_ve_scores is not None and len(user_ve_scores) == 0:
            user_ve_scores = None
        if variant_effect_sources is not None and len(variant_effect_sources) == 0:
            variant_effect_sources = None

        # Get the full universe of variants for query criteria. The universe
        # is limited to those for which we have labels.
        variant_universe_pks_df = self._variant_effect_label_repo.get(
            task_code, variant_query_criteria
        )[VARIANT_PK_COLUMNS]

        # if user has specified variant scores then restrict user
        # variants to those in the universe(i.e. those for which we have
        # labels) and then reset the universe to the filtered user variant list
        if user_ve_scores is not None:
            if column_name_map is not None and len(column_name_map) > 0:
                user_ve_scores = user_ve_scores.rename(columns=column_name_map)
            self.validate_user_ve_scores(user_ve_scores)
            user_ve_scores = user_ve_scores.merge(
                variant_universe_pks_df, how="inner", on=VARIANT_PK_COLUMNS
            )
            user_ve_scores = user_ve_scores.assign(SCORE_SOURCE=user_vep_name)
            variant_universe_pks_df = user_ve_scores[VARIANT_PK_COLUMNS]

        if (
            user_ve_scores is not None
            and not include_variant_effect_sources
            and variant_effect_sources is None
        ):
            # We only do the analysis against the user scores.
            # We don't use any system veps.
            analysis_ve_scores_df = user_ve_scores[VARIANT_EFFECT_SCORE_COLS]
            analysis_labels_df = self._variant_effect_label_repo.get(
                task_code,
                VEQueryCriteria(variant_ids=user_ve_scores[VARIANT_PK_COLUMNS]),
            )
        else:
            # We include system veps in the analysis
            # get the system vep scores based upon selection critera
            system_ve_scores_df = self._variant_effect_score_repo.get(
                task_code,
                variant_effect_sources,
                include_variant_effect_sources,
                variant_query_criteria,
            )

            # First restrict the set of system vep scores to the variants
            # in the universe. Then compute how many variants there are for
            # each vep. Then we only keep the vep scores for veps where
            # the variant count is above the vep_min_overlap_count
            vep_min_overlap_count = (
                len(variant_universe_pks_df) * vep_min_overlap_percent * 0.01
            )
            system_ve_scores_df = filter_dataframe_by_list(
                system_ve_scores_df, variant_universe_pks_df, VARIANT_PK_COLUMNS
            )
            retained_veps = []
            system_ve_scores_count_by_vep = system_ve_scores_df.groupby(
                "SCORE_SOURCE"
            ).size()
            retained_veps = system_ve_scores_count_by_vep[
                system_ve_scores_count_by_vep >= vep_min_overlap_count
            ].index
            system_ve_scores_df = filter_dataframe_by_list(
                system_ve_scores_df, retained_veps, "SCORE_SOURCE"
            )

            # Now for each variant we compute how many veps for which we have
            # scores. We then retain only those variants where the number of
            # veps is above variant_vep_retention_count
            variant_vep_retention_count = (
                len(retained_veps) * variant_vep_retention_percent * 0.01
            )

            system_ve_scores_count_by_var = system_ve_scores_df.groupby(
                VARIANT_PK_COLUMNS
            ).size()
            retained_variants = system_ve_scores_count_by_var[
                system_ve_scores_count_by_var >= variant_vep_retention_count
            ].reset_index()

            # If user specified scores append them to the system ones. Then
            # merge in the labels for all of the variants.
            if user_ve_scores is not None:
                analysis_ve_scores_df = pd.concat(
                    [
                        system_ve_scores_df[VARIANT_EFFECT_SCORE_COLS],
                        user_ve_scores[VARIANT_EFFECT_SCORE_COLS],
                    ]
                )
            else:
                analysis_ve_scores_df = system_ve_scores_df[VARIANT_EFFECT_SCORE_COLS]
            analysis_ve_scores_df = filter_dataframe_by_list(
                analysis_ve_scores_df, retained_variants, VARIANT_PK_COLUMNS
            )
            analysis_labels_df = self._variant_effect_label_repo.get(
                task_code, VEQueryCriteria(variant_ids=retained_variants)
            )
        analysis_ve_scores_labels_df = analysis_ve_scores_df.merge(
            analysis_labels_df, how="inner", on=VARIANT_PK_COLUMNS
        )
        return analysis_ve_scores_labels_df

    def _compute_pr(
        self,
        agg_level: str,
        grouped_ve_scores_labels,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:

        score_sources = []
        aucs = []
        genes = []
        precision_recall_score_sources = []
        precision_recall_genes = []
        precisions = []
        recalls = []
        thresholds_list = []
        include_gene_metrics = agg_level == "vep_gene"
        for score_source_gene, scores_labels_df in grouped_ve_scores_labels:
            score_source = score_source_gene[0]
            score_sources.append(score_source)
            if include_gene_metrics:
                genes.append(score_source_gene[1])
            prs, recs, thresholds = precision_recall_curve(
                scores_labels_df["BINARY_LABEL"], scores_labels_df["RANK_SCORE"]
            )
            aucs.append(auc(recs, prs))
            precision_recall_score_sources.extend([score_source] *
                                                  len(thresholds))
            if include_gene_metrics:
                precision_recall_genes.extend([score_source_gene[1]] *
                                              len(thresholds))
            precisions.extend(prs[:-1])
            recalls.extend(recs[:-1])
            thresholds_list.extend(thresholds)
        pr_df = pd.DataFrame(
            {"SCORE_SOURCE": score_sources, "PR_AUC": aucs}
        )
        pr_curve_coords_df = pd.DataFrame(
            {
                "SCORE_SOURCE": precision_recall_score_sources,
                "PRECISION": precisions,
                "RECALL": recalls,
                "THRESHOLD": thresholds_list,
            }
        )
        if include_gene_metrics:
            pr_df["GENE_SYMBOL"] = genes
            pr_curve_coords_df["GENE_SYMBOL"] = precision_recall_genes
        return pr_df, pr_curve_coords_df

    def _compute_roc(
        self,
        agg_level: str,
        grouped_ve_scores_labels,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:

        score_sources = []
        aucs = []
        genes = []
        fpr_tpr_score_sources = []
        fpr_tpr_genes = []
        false_positive_rates = []
        true_positive_rates = []
        thresholds_list = []
        exceps = []
        include_gene_metrics = agg_level == "vep_gene"
        for score_source_gene, scores_labels_df in grouped_ve_scores_labels:
            try:
                score_source = score_source_gene[0]
                score_sources.append(score_source)
                if include_gene_metrics:
                    genes.append(score_source_gene[1])
                if (
                    len(scores_labels_df) == sum(
                        scores_labels_df["BINARY_LABEL"])
                    or sum(scores_labels_df["BINARY_LABEL"]) == 0
                ):
                    raise Exception(
                        "Cannot compute roc metrics because all "
                        "labels have same value"
                    )
                auc = roc_auc_score(
                    scores_labels_df["BINARY_LABEL"], scores_labels_df[
                        "RANK_SCORE"]
                )
                fpr, tpr, thresholds = roc_curve(
                    scores_labels_df["BINARY_LABEL"], scores_labels_df[
                        "RANK_SCORE"]
                )
                aucs.append(auc)
                fpr_tpr_score_sources.extend([score_source] * len(fpr))
                if include_gene_metrics:
                    fpr_tpr_genes.extend([score_source_gene[1]] * len(fpr))
                false_positive_rates.extend(fpr)
                true_positive_rates.extend(tpr)
                thresholds_list.extend(thresholds)
                exceps.append(np.nan)
            except Exception as e:
                aucs.append(np.nan)
                exceps.append(str(e))

        roc_df = pd.DataFrame(
            {
                "SCORE_SOURCE": score_sources,
                "ROC_AUC": aucs,
                "EXCEPTION": exceps,
            }
        )
        roc_curve_coords_df = pd.DataFrame(
            {
                "SCORE_SOURCE": fpr_tpr_score_sources,
                "FALSE_POSITIVE_RATE": false_positive_rates,
                "TRUE_POSITIVE_RATE": true_positive_rates,
                "THRESHOLD": thresholds_list,
            }
        )
        if include_gene_metrics:
            roc_df["GENE_SYMBOL"] = genes
            roc_curve_coords_df["GENE_SYMBOL"] = fpr_tpr_genes
        return roc_df, roc_curve_coords_df

    @staticmethod
    def _compute_general_metrics(group) -> pd.Series:
        return pd.Series(
            {
                "NUM_VARIANTS": len(group),
                "NUM_POSITIVE_LABELS": group["BINARY_LABEL"].sum(),
                "NUM_NEGATIVE_LABELS": (group["BINARY_LABEL"] ^ 1).sum(),
            }
        )

    @staticmethod
    def _compute_num_unique_variants(group) -> pd.Series:
        return pd.Series(
                {
                    "NUM_UNIQUE_VARIANTS": group[VARIANT_PK_COLUMNS].
                    drop_duplicates().shape[0]
                })

    def _add_info_to_metric_dataframes(self, *dfs):  # -> list(pd.DataFrame):
        return_dfs = []
        for df in dfs:
            if df is not None:
                df = df.merge(
                    self._variant_effect_source_repo.get_all()[["CODE", "NAME"]],
                    how="left",
                    left_on="SCORE_SOURCE",
                    right_on="CODE",
                )
                df.loc[df["NAME"].isna(), "NAME"] = df.loc[
                    df["NAME"].isna(), "SCORE_SOURCE"
                ]
                df.rename(columns={"NAME": "SOURCE_NAME"}, inplace=True)
                df.drop(columns="CODE", inplace=True)
            return_dfs.append(df)
        return return_dfs

    def _compute_mwu(
        self,
        agg_level: str,
        grouped_ve_scores_labels,
    ) -> pd.DataFrame:

        score_sources = []
        neg_log10_mwu_pvals = []
        genes = []
        exceps = []
        include_gene_metrics = agg_level == "vep_gene"
        for score_source_gene, scores_labels_df in grouped_ve_scores_labels:
            try:
                score_source = score_source_gene[0]
                score_sources.append(score_source)
                if include_gene_metrics:
                    genes.append(score_source_gene[1])
                positive_scores = np.array(
                    scores_labels_df.query("BINARY_LABEL == 1")["RANK_SCORE"]
                )
                negative_scores = np.array(
                    scores_labels_df.query("BINARY_LABEL == 0")["RANK_SCORE"]
                )
                if 0 in [len(positive_scores), len(negative_scores)]:
                    raise Exception(
                        "Cannot compute mann-whitney u values "
                        + "because all labels have same value"
                    )
                neg_log10_mwu_pval = -np.log10(
                    stats.mannwhitneyu(negative_scores, positive_scores).pvalue
                )
                neg_log10_mwu_pvals.append(neg_log10_mwu_pval)
                exceps.append(np.nan)
            except Exception as e:
                neg_log10_mwu_pvals.append(np.nan)
                exceps.append(str(e))

        mwu_df = pd.DataFrame(
            {
                "SCORE_SOURCE": score_sources,
                "NEG_LOG10_MWU_PVAL": neg_log10_mwu_pvals,
                "EXCEPTION": exceps,
            }
        )
        if include_gene_metrics:
            mwu_df["GENE_SYMBOL"] = genes
        return mwu_df

    def _compute_metrics(
        self,
        task_code: str,
        ve_scores_labels_df: pd.DataFrame,
        compute_gene_metrics: bool,
        metrics: list[str],
        list_variants: bool = False,
    ):

        grouped_ve_scores_labels = ve_scores_labels_df.groupby(
            ["SCORE_SOURCE"])
        general_metrics_df = grouped_ve_scores_labels.apply(
            self._compute_general_metrics, include_groups=False
        ).reset_index()
        gene_general_metrics_df = None
        if compute_gene_metrics:
            grouped_ve_gene_scores_labels = ve_scores_labels_df.groupby(
                ["SCORE_SOURCE", "GENE_SYMBOL"]
            )
            gene_general_metrics_df = grouped_ve_gene_scores_labels.apply(
                self._compute_general_metrics, include_groups=False
            ).reset_index()
        roc_df = None
        roc_curve_coords_df = None
        pr_df = None
        pr_curve_coords_df = None
        mwu_df = None
        gene_roc_df = None
        gene_roc_curve_coords_df = None
        gene_pr_df = None
        gene_pr_curve_coords_df = None
        gene_mwu_df = None
        if "roc" in metrics:
            roc_df, roc_curve_coords_df = self._compute_roc(
                "vep", grouped_ve_scores_labels)
            if compute_gene_metrics:
                gene_roc_df, gene_roc_curve_coords_df = self._compute_roc(
                    "vep_gene", grouped_ve_gene_scores_labels)
        if "pr" in metrics:
            pr_df, pr_curve_coords_df = self._compute_pr(
                "vep",
                grouped_ve_scores_labels)
            if compute_gene_metrics:
                gene_pr_df, gene_pr_curve_coords_df = self._compute_pr(
                    "vep_gene", grouped_ve_gene_scores_labels)
        if "mwu" in metrics:
            mwu_df = self._compute_mwu("vep", grouped_ve_scores_labels)
            if compute_gene_metrics:
                gene_mwu_df = self._compute_mwu(
                    "vep_gene", grouped_ve_gene_scores_labels)
        if list_variants:
            included_variants_df = ve_scores_labels_df[
                ["SCORE_SOURCE"] + VARIANT_PK_COLUMNS
            ]
        else:
            included_variants_df = None
        metric_dataframes = self._add_info_to_metric_dataframes(
            general_metrics_df, roc_df, pr_df, mwu_df, gene_general_metrics_df,
            gene_roc_df, gene_pr_df,
            gene_mwu_df
        )
        metric_dataframes.extend(
            [roc_curve_coords_df, pr_curve_coords_df, gene_roc_curve_coords_df,
             gene_pr_curve_coords_df, included_variants_df]
        )
        return metric_dataframes

    def _filter_scores_and_labels_by_top_genes(
        self, scores_and_labels_df: pd.DataFrame, num_top_genes: int
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        gene_unique_variant_counts_df = scores_and_labels_df.groupby(
            ["GENE_SYMBOL"]).apply(
                self._compute_num_unique_variants, include_groups=False
        ).reset_index()
        gene_unique_variant_counts_df = (
            gene_unique_variant_counts_df.sort_values(
                by="NUM_UNIQUE_VARIANTS", ascending=False
            ).iloc[:num_top_genes]
        )
        scores_and_labels_df = scores_and_labels_df.merge(
            gene_unique_variant_counts_df, how="inner", on="GENE_SYMBOL"
        )
        return gene_unique_variant_counts_df, scores_and_labels_df

    def compute_metrics(
        self,
        task_code: str,
        user_ve_scores: pd.DataFrame = None,
        user_vep_name: str = "USER",
        compute_gene_metrics: bool = False,
        num_top_genes: int = None,
        column_name_map: dict = None,
        variant_effect_sources: list[str] = None,
        include_variant_effect_sources: bool = True,
        variant_query_criteria: VEQueryCriteria = None,
        vep_min_overlap_percent: float = 0,
        variant_vep_retention_percent: float = 0,
        metrics: str | list[str] = ["roc", "pr", "mwu"],
        list_variants: bool = False,
    ) -> VEAnalysisResult:
        """
        Generates performance metrics for an optional user supplied set of
        vep scores and for system supplied vep's. If the user doesn't provide
        vep scores, will only generate metrics for system veps. Returns
        an object containing all the metrics which can then be used to
        generate plots, reports, or csv data files.

        Parameters
        ----------
        task_code : str
            Task code
        user_ve_scores : DataFrame, optional
            An optional dataframe of user variant effect prediction
            scores. Expected to have the following columns:
            GENOME_ASSEMBLY, CHROMOSOME, POSITION,
            REFERENCE_NUCLEOTIDE, ALTERNATE_NUCLEOTIDE, RANK_SCORE.
            The GENOME_ASSEMBLY must be hg38 in current release.
            RANK_SCORE is a numeric prediction score. It does not have
            to be standardized or normalized.
        user_vep_name : str, optional
            If user_ve_scores are provided, then this is the label to
            be used for them in the analysis output.
        compute_gene_metrics: bool
            Whether to compute vep/gene level metrics along with vep level
            metrics.
        num_top_genes: int
            If specified, only consider the top N genes by number of
            unique variants that satisfy the selection criteria indicated
            by the other parameters.
        agg_level : str
            The level at which to aggregate the metrics. Can be either
            "vep" or "vep_gene". If "vep", then the metrics are aggregated
            at the vep level. If "vep_gene", then the metrics are aggregated
            at the vep/gene level.
        column_name_map : dict, optional
            If the column names in user_ve_scores are not the expected
            names, then this maps the column names to the expected names.
        variant_effect_sources : list, optional
            If specified it would restrict the analysis to the
            system supplied vep's in this list.
        include_variant_effect_sources : bool, optional
            If variant_effect_sources is specified, indicates whether to
            limit the analysis to the system supplied vep's specified or to
            exclude the system supplied vep's specified.
            If variant_effect_sources is not specified a value of False
            indicates that all system variant effect sources should
            be excluded from the analysis.
        variant_query_criteria : VEQueryCriteria, optional
            See description of VEQueryCriteria in model package.
            Specifies criteria that would limit the set of variants
            to be included in the analysis.
        vep_min_overlap_percent : float
            In order for a system supplied vep to be included in the
            analysis the set of variants for which it has scores must
            overlap the variants in user_ve_scores by at least this
            percentage amount. If the user_ve_scores is not specified,
            then it must overlap the entire universe of variants in
            the system for which we have labels by at least this
            percentage amount. A value of 0 means that it can overlap
            by any amount to be included.
        variant_vep_retention_percent : float
            In order for a variant to be included in the analysis
            there must exist scores for the variant in at least
            this percent of the system vep's included in the analysis
            based on the value of vep_min_overlap_percent. For example,
            if a value of 50 is specified and 10 system vep's qualified
            for inclusion, then a variant must be in at least 5 veps
            to be included in the analysis.
        metrics : str or list[str]
            Specifies which metrics to compute. Can be a string
            indicating a single metric or a list of strings for
            multiple metrics. The metrics are: roc, pr, mwu.
        list_variants: bool
            Include the list of variants
            that were included in the analysis in the return result
            object. There is a separate list for the user variants
            as well as for each system vep.

        Returns
        -------
        VEAnalysisResult
            Object containing computed metrics
        """

        validate_query_criteria(variant_query_criteria)
        if not compute_gene_metrics:
            num_top_genes = None

        scores_and_labels_df = self.get_analysis_scores_and_labels(
            task_code,
            user_ve_scores,
            user_vep_name,
            column_name_map,
            variant_effect_sources,
            include_variant_effect_sources,
            variant_query_criteria,
            vep_min_overlap_percent,
            variant_vep_retention_percent,
        )

        if num_top_genes is not None:
            (
                gene_unique_variant_counts_df,
                scores_and_labels_df
            ) = self._filter_scores_and_labels_by_top_genes(
                scores_and_labels_df, num_top_genes
            )
        else:
            gene_unique_variant_counts_df = None

        if type(metrics) is str:
            metrics = [metrics]
        (
            general_metrics_df,
            roc_df,
            pr_df,
            mwu_df,
            gene_general_metrics_df,
            gene_roc_df,
            gene_pr_df,
            gene_mwu_df,
            roc_curve_coords_df,
            pr_curve_coords_df,
            gene_roc_curve_coords_df,
            gene_pr_curve_coords_df,
            included_variants_df,
        ) = self._compute_metrics(
            task_code, scores_and_labels_df, compute_gene_metrics, metrics,
            list_variants
        )
        num_variants = len(
            scores_and_labels_df[VARIANT_PK_COLUMNS].drop_duplicates())
        num_user_variants = None if user_ve_scores is None else len(
            user_ve_scores)
        return VEAnalysisResult(
            num_variants,
            num_user_variants,
            user_vep_name,
            general_metrics_df,
            roc_df,
            pr_df,
            mwu_df,
            gene_general_metrics_df,
            gene_roc_df,
            gene_pr_df,
            gene_mwu_df,
            roc_curve_coords_df,
            pr_curve_coords_df,
            gene_roc_curve_coords_df,
            gene_pr_curve_coords_df,
            included_variants_df,
            gene_unique_variant_counts_df
        )
