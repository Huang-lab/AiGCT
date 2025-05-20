"""
Classes that represent model objects.

Model objects are containers for data. They generally do not have
behavior associated with them.
"""

from dataclasses import dataclass
import pandas as pd
from typing import List, Dict


@dataclass
class VariantId:
    """
    Model object that represents a variant id.

    Attributes
    ----------
    genome_assembly : str
        genome assembly symbol, i.e. hg38
    chromosome : str
        chromosome
    position : int
        position
    reference_nucleotide : str
        reference nucleotide
    alternate_nucleotide : str
        alternate_nucleotide
    """

    genome_assembly: str
    chromosome: str
    position: int
    reference_nucleotide: str
    alternate_nucleotide: str


@dataclass
class VariantEffectSource:
    """
    Model object that represents a variant effect source.

    Attributes
    ----------
    code : str
        A unique code that identifies the variant effect source
    name : str
        A unique name of the source
    source_type : str
        i.e. VEP
    description : str
        Description
    """

    code: str
    name: str
    source_type: str
    description: str


@dataclass
class VariantFilter:
    """
    Model object that represents named variant filter query. The
    filter criteria consists either of a list of genes or a list
    of variant id's or both.

    Attributes
    ----------
    filter : Series
        A series with the the unique code, name, description of the
        filter.
    filter_genes : DataFrame
        A dataframe of gene symbols associated with the filter.
        If None then there filter_variants must not be None.
    filter_variants : DataFrame
        A dataframe of variant id's associated with the filter.
        If None then there filter_genes must not be None.
    """

    filter: pd.Series
    filter_genes: pd.DataFrame
    filter_variants: pd.DataFrame


@dataclass
class VEQueryCriteria:
    """
    Model object that represents variant query criteria.

    Attributes
    ----------
    gene_symbols : list or DataFrame, optional
        List of gene symbols
    include_genes : bool, optional
        If gene_symbols is provided, indicates whether to limit variants
        to associated with those gene_symbols to exclude variants
        associated with the gene_symbols.
    variant_ids : DataFrame, optional
        List of variant ids. The dataframe is expected to have the
        following columns:
        GENOME_ASSEMBLY, CHROMOSOME, POSITION,
        REFERENCE_NUCLEOTIDE, ALTERNATE_NUCLEOTIDE
        If the column names are different specify a value for
        column_name_map mapping the column names to the expected names.
    include_variant_ids : bool, optional
        If variant_ids is provided, indicates whether to limit variants
        to the variant_ids provided or to fetch all variants but those
        in variant_ids
    column_name_map : Dict, optional
        A dictionary that maps the column names in variant_ids to the
        expected column names.
    allele_frequency_operator : str, optional
        If allele_frequency is provided, this is one of "eq", "gt",
        "lt", "ge", "le". i.e. limit variants to those whose
        allele_frequency is equal to, greater than, etc. the
        allele_frequency.
    allele_frequency : float, optional
        Used in conjunction to allele_frequency_operator to limit variants
        to those meeting a certain allele_frequency criteria.
    filter_names : str | list[str], optional
        The name(s) of a system filter that can be used to limit the variants
        returned. If more than one is given then the filters are combined
        using a logical OR.
    """

    gene_symbols: list[str] | pd.DataFrame | pd.Series = None
    include_genes: bool = True
    variant_ids: pd.DataFrame = None
    include_variant_ids: bool = True
    column_name_map: Dict = None
    allele_frequency_operator: str = "="
    allele_frequency: float = None
    filter_names: str | list[str] = None


@dataclass
class VEAnalysisResult:
    """
    Represents the result of calling VEAnalyzer.compute_metrics.

    Attributes
    ----------
    num_variants_included : int
        The total number of unique variants included in the analysis
        across all veps.
    num_user_variants : int
        The number of user supplied variants included in the analysis
    user_vep_name : str
        Name of user vep
    general_metrics : DataFrame
        Has the following columns:
        SCORE_SOURCE - Short unique vep identifier
        NUM_VARIANTS, NUM_POSITIVE_LABELS, NUM_NEGATIVE_LABELS,
        SOURCE_NAME - Name of vep
    roc_metrics : DataFrame, optional
        Roc metrics with columns: SCORE_SOURCE,
        ROC_AUC, EXCEPTION, SOURCE_NAME
        EXCEPTION would store an exception message in the event the
        roc could not be computed for that vep.
    pr_metrics : DataFrame, optional
        Precision/Recall metrics containing columns: SCORE_SOURCE,
        PR_AUC, SOURCE_NAME
    mwu_metrics : DataFrame, optional
        Mann-Whitney U metrics containing columns: SCORE_SOURCE,
        NEG_LOG10_MWU_PVAL, SOURCE_NAME
    gene_general_metrics : DataFrame, optional
        Gene-level general metrics with columns: SCORE_SOURCE,
        GENE_SYMBOL, NUM_VARIANTS, NUM_POSITIVE_LABELS,
        NUM_NEGATIVE_LABELS, SOURCE_NAME
    gene_roc_metrics : DataFrame, optional
        Gene-level ROC metrics with columns: SCORE_SOURCE,
        GENE_SYMBOL, ROC_AUC, EXCEPTION, SOURCE_NAME
    gene_pr_metrics : DataFrame, optional
        Gene-level precision/recall metrics with columns: SCORE_SOURCE,
        GENE_SYMBOL, PR_AUC, SOURCE_NAME
    gene_mwu_metrics : DataFrame, optional
        Gene-level Mann-Whitney U metrics with columns: SCORE_SOURCE,
        GENE_SYMBOL, NEG_LOG10_MWU_PVAL, SOURCE_NAME
    roc_curve_coordinates : DataFrame, optional
        Columns: SCORE_SOURCE,
        FALSE_POSITIVE_RATE, TRUE_POSITIVE_RATE, THRESHOLD
    pr_curve_coordinates : DataFrame, optional
        Columns: SCORE_SOURCE,
        PRECISION, RECALL, THRESHOLD
    gene_roc_curve_coordinates : DataFrame, optional
        Gene-level ROC curve coordinates with columns: SCORE_SOURCE,
        GENE_SYMBOL, FALSE_POSITIVE_RATE, TRUE_POSITIVE_RATE, THRESHOLD
    gene_pr_curve_coordinates : DataFrame, optional
        Gene-level precision/recall curve coordinates with columns:
        SCORE_SOURCE, GENE_SYMBOL, PRECISION, RECALL, THRESHOLD
    variants_included : DataFrame, optional
        List of variants included for each vep included the user vep.
        Columns:
        SCORE_SOURCE, GENOME_ASSEMBLY, CHROMOSOME, POSITION,
        REFERENCE_NUCLEOTIDE, ALTERNATE_NUCLEOTIDE
    gene_unique_variant_counts_df : DataFrame, optional
        Count of unique variants per gene across all vepswith columns:
        GENE_SYMBOL, NUM_UNIQUE_VARIANTS
    """

    num_variants_included: int
    num_user_variants: int
    user_vep_name: str
    general_metrics: pd.DataFrame
    roc_metrics: pd.DataFrame
    pr_metrics: pd.DataFrame
    mwu_metrics: pd.DataFrame
    gene_general_metrics: pd.DataFrame
    gene_roc_metrics: pd.DataFrame
    gene_pr_metrics: pd.DataFrame
    gene_mwu_metrics: pd.DataFrame
    roc_curve_coordinates: pd.DataFrame
    pr_curve_coordinates: pd.DataFrame
    gene_roc_curve_coordinates: pd.DataFrame
    gene_pr_curve_coordinates: pd.DataFrame
    variants_included: pd.DataFrame
    gene_unique_variant_counts_df: pd.DataFrame


@dataclass
class TaskPkViolations:
    dups_found: bool
    variant_effect_label_dups: pd.DataFrame
    variant_effect_score_dups: pd.DataFrame
    variant_filter_dups: pd.DataFrame
    variant_filter_gene_dups: pd.DataFrame
    variant_filter_variant_dups: pd.DataFrame


@dataclass
class PkViolations:
    dups_found: bool
    variant_dups: pd.DataFrame
    variant_effect_source_dups: pd.DataFrame
    task_violations: dict[str, TaskPkViolations]
