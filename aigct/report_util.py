import pandas as pd


class GeneMetricSorter:
    """
    Sort the gene metrics dataframe by variant effect source ascending
    and descending by the number of unique variants
    per gene. This will ensure that the gene metrics are presented
    in a consistent order across all gene metric dataframes.
    """

    def __init__(self, gene_variant_counts: pd.DataFrame,
                 num_top_genes: int = None):
        """
        Parameters
        ----------
        gene_variant_counts : pd.DataFrame
            Dataframe with columns GENE_SYMBOL, NUM_UNIQUE_VARIANTS
        num_top_genes : int, optional
            If specified, only consider the top N genes by number of
            unique variants.
        """
        if num_top_genes is not None:
            self._gene_variant_counts = gene_variant_counts.sort_values(
                by="NUM_UNIQUE_VARIANTS", ascending=False   
            ).iloc[:num_top_genes]
        else:
            self._gene_variant_counts = gene_variant_counts

    def sort_gene_metrics(self, gene_metrics: pd.DataFrame):
        cols = gene_metrics.columns
        if "SCORE_SOURCE" in cols:
            score_source_col = "SCORE_SOURCE"
        else:
            score_source_col = "SOURCE_NAME"
        return gene_metrics.merge(self._gene_variant_counts, how="inner",
                                  on="GENE_SYMBOL").sort_values(
                                     by=[score_source_col,
                                         "NUM_UNIQUE_VARIANTS"],
                                     ascending=[True, False]
                                     )[cols]
        

