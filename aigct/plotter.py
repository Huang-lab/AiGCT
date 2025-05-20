import matplotlib.pyplot as plt
# from IPython.display import HTML
from IPython.display import display
# import seaborn as sns
# import plotly.express as px
import pandas as pd
import os
from .model import (
    VEAnalysisResult
)
from .util import Config
from .date_util import now_str_compact
from .file_util import (
    unique_file_name,
    create_folder
)
from .plot_util import (barchart, create_feature_palette,
                        get_colors, heatmap)
from .report_util import GeneMetricSorter


class VEAnalysisPlotter:
    """Plot results of an analysis"""

    def __init__(self, config: Config):
        self._config = config
        self._roc_pr_config = config.roc_pr_line
        self._mwu_config = config.mwu_bar

    def _plot_roc_curves(self, aucs: pd.DataFrame, user_vep_name: str,
                         curve_coords: pd.DataFrame,
                         batch_no: int, num_batches: int,
                         ves_color_palette: dict,
                         file_name: str = None):
        title = self._roc_pr_config.roc_title + (
            f" ({batch_no+1} of {num_batches})" if num_batches > 1 else "")
        plt.figure(figsize=(10, 8))
        for i, ve_auc in aucs.iterrows():
            if pd.isna(ve_auc['ROC_AUC']):
                continue
            ve_curve_coords = curve_coords[curve_coords['SCORE_SOURCE'] ==
                                           ve_auc['SCORE_SOURCE']].sort_values(
                'THRESHOLD', ascending=False)
            label = ve_auc['SOURCE_NAME'] + \
                ' (AUC=' + str(round(ve_auc['ROC_AUC'], 4)) + ')'
            plt.plot(
                ve_curve_coords['FALSE_POSITIVE_RATE'],
                ve_curve_coords['TRUE_POSITIVE_RATE'],
                label=label,
                color=ves_color_palette[ve_auc['SOURCE_NAME']],
                lw=self._config.line_width,
                linestyle=self._config.line_style)

        plt.xlabel('False Positive Rate',
                   fontsize=self._roc_pr_config.x_axis_font_size)
        plt.ylabel('True Positive Rate',
                   fontsize=self._roc_pr_config.y_axis_font_size)
        plt.tick_params(axis='both', labelsize=self._roc_pr_config.label_size)
        plt.title(title, fontsize=self._roc_pr_config.title_font_size)
        legend = plt.legend(loc="lower right",
                            fontsize=self._roc_pr_config.legend_font_size)
        for line in legend.get_lines():
            line.set_linewidth(self._roc_pr_config.legend_line_width)
        if file_name is None:
            plt.show()
            return
        plt.savefig(file_name, dpi=self._config.file_dpi,
                    format='png',
                    bbox_inches=self._config.bbox_inches)
        plt.savefig(file_name.replace(".png", ".svg"),
                    format='svg',
                    bbox_inches=self._config.bbox_inches)

    def _plot_pr_curves(self, aucs: pd.DataFrame, user_vep_name: str,
                        curve_coords: pd.DataFrame,
                        batch_no: int, num_batches: int,
                        ves_color_palette: dict, file_name: str = None):
        title = self._roc_pr_config.pr_title + \
            (f" ({batch_no+1} of {num_batches})" if num_batches > 1 else "")
        plt.figure(figsize=(10, 8))
        for i, ve_auc in aucs.iterrows():
            ve_curve_coords = curve_coords[curve_coords['SCORE_SOURCE'] ==
                                           ve_auc['SCORE_SOURCE']].sort_values(
                'THRESHOLD', ascending=False)
            label = ve_auc['SOURCE_NAME'] + \
                ' (AUC=' + str(round(ve_auc['PR_AUC'], 4)) + ')'
            plt.plot(
                ve_curve_coords['RECALL'][1:],
                ve_curve_coords['PRECISION'][1:],
                label=label,
                color=ves_color_palette[ve_auc['SOURCE_NAME']],
                lw=self._config.line_width,
                linestyle=self._config.line_style)

        plt.xlabel('Recall', fontsize=self._roc_pr_config.x_axis_font_size)
        plt.ylabel('Precision', fontsize=self._roc_pr_config.y_axis_font_size)
        plt.tick_params(axis='both', labelsize=self._roc_pr_config.label_size)
        plt.title(title, fontsize=self._roc_pr_config.title_font_size)
        plt.ylim([0, 1.05])
        legend = plt.legend(loc="lower center",
                            fontsize=self._roc_pr_config.legend_font_size)
        for line in legend.get_lines():
            line.set_linewidth(self._roc_pr_config.legend_line_width)
        if file_name is None:
            plt.show()
            return
        plt.savefig(file_name, dpi=self._config.file_dpi,
                    format='png',
                    bbox_inches=self._config.bbox_inches)
        plt.savefig(file_name.replace(".png", ".svg"),
                    format='svg',
                    bbox_inches=self._config.bbox_inches)

    def _display_mwu_table(self, results: VEAnalysisResult,
                           file_name: str = None):
        table_output = results.general_metrics.merge(
            results.mwu_metrics, how="inner", suffixes=(None, "_y"),
            on='SCORE_SOURCE')[[
                'SOURCE_NAME', 'NUM_VARIANTS', 'NUM_POSITIVE_LABELS',
                'NUM_NEGATIVE_LABELS', 'NEG_LOG10_MWU_PVAL', 'EXCEPTION']]
        table_output = table_output.sort_values('NEG_LOG10_MWU_PVAL',
                                                ascending=False)
        style = table_output.style.set_properties(
            subset=["NUM_VARIANTS",
                    "NUM_POSITIVE_LABELS", "NUM_NEGATIVE_LABELS",
                    "NEG_LOG10_MWU_PVAL"], **{"text-align": "right"})
        style = style.hide().relabel_index(
             ['VEP', "Variant Total", "Positive Labels", "Negative Labels",
              "MWU -log10(pval)", "Exception"],
             axis=1).set_caption("Mann-Whitney U -log10(p value)")
        if file_name:
            style.to_html(file_name)
        else:
            display(style)

    def _display_pr_table(self, results: VEAnalysisResult,
                          file_name: str = None):
        table_output = results.general_metrics.merge(
            results.pr_metrics, how="inner", suffixes=(None, "_y"),
            on='SCORE_SOURCE')[[
                'SOURCE_NAME', 'NUM_VARIANTS', 'NUM_POSITIVE_LABELS',
                'NUM_NEGATIVE_LABELS', 'PR_AUC']]
        table_output = table_output.sort_values('PR_AUC', ascending=False)
        style = table_output.style.set_properties(
            subset=["NUM_VARIANTS",
                    "NUM_POSITIVE_LABELS", "NUM_NEGATIVE_LABELS",
                    "PR_AUC"], **{"text-align": "right"})
        style = style.hide().relabel_index(
             ['VEP', "Variant Total", "Positive Labels", "Negative Labels",
              "PR AUC"],
             axis=1).set_caption("Precision/Recall")
        style = style.set_properties(
            subset=["NUM_VARIANTS",
                    "NUM_POSITIVE_LABELS", "NUM_NEGATIVE_LABELS",
                    "PR_AUC"], **{"text-align": "right"})
        if file_name:
            style.to_html(file_name)
        else:
            display(style)

    def _display_roc_table(self, results: VEAnalysisResult,
                           file_name: str = None):
        table_output = results.general_metrics.merge(
            results.roc_metrics, how="inner", suffixes=(None, "_y"),
            on='SCORE_SOURCE')[[
                'SOURCE_NAME', 'NUM_VARIANTS', 'NUM_POSITIVE_LABELS',
                'NUM_NEGATIVE_LABELS', 'ROC_AUC']]
        table_output = table_output.sort_values('ROC_AUC', ascending=False)
        style = table_output.style.set_properties(
            subset=["NUM_VARIANTS",
                    "NUM_POSITIVE_LABELS", "NUM_NEGATIVE_LABELS",
                    "ROC_AUC"], **{"text-align": "right"})
        style = style.hide().relabel_index(
             ['VEP', "Variant Total", "Positive Labels", "Negative Labels",
              "ROC AUC"],
             axis=1).set_caption("ROC")
        if file_name:
            style.to_html(file_name)
        else:
            display(style)

    def _plot_mwu_bar(self, mwus: pd.DataFrame,
                      batch_no: int, num_batches: int,
                      file_name: str = None, palette=None):
        title = self._mwu_config.title + \
            (f" ({batch_no+1} of {num_batches})" if num_batches > 1 else "")
        config = self._mwu_config
        barchart(mwus, 'SOURCE_NAME', 'NEG_LOG10_MWU_PVAL',
                 palette=palette,
                 y_label='Mann-Whitney U log10(p value)',
                 filename=file_name, title=title,
                 title_fontsize=config.title_font_size,
                 y_label_fontsize=config.y_label_font_size,
                 x_label_fontsize=config.x_label_font_size,
                 label_size=config.label_size,
                 xtick_rotation=config.xtick_rotation,
                 xtick_rotation_mode=config.xtick_rotation_mode,
                 file_dpi=self._config.file_dpi,
                 bbox_inches=self._config.bbox_inches)

    def plot_pr_results(self, results: VEAnalysisResult,
                        ves_color_palette: dict,
                        dir: str = None):
        num_curves_per_plot = self._roc_pr_config.num_curves_per_plot
        plot_batches = []
        pr_metrics = results.pr_metrics.sort_values('PR_AUC',
                                                    ascending=False)
        for idx in range(0, len(pr_metrics), num_curves_per_plot):
            batch = pr_metrics.iloc[idx:idx+num_curves_per_plot]
            plot_batches.append(batch)
        pr_curves_file_name = None if dir is None else os.path.join(
            dir, "pr_curves_")
        num_batches = len(plot_batches)
        for batch_no, batch in enumerate(plot_batches):
            batch_file_name = None if pr_curves_file_name is None else \
                pr_curves_file_name + str(batch_no) + ".png"
            self._plot_pr_curves(batch, results.user_vep_name,
                                 results.pr_curve_coordinates,
                                 batch_no, num_batches,
                                 ves_color_palette, batch_file_name)
        pr_table_file_name = None if dir is None else os.path.join(
            dir, "pr_table" + ".html")
        self._display_pr_table(results, pr_table_file_name)

    def plot_roc_results(self, results: VEAnalysisResult,
                         ves_color_palette: dict,
                         dir: str = None):
        num_curves_per_plot = self._roc_pr_config.num_curves_per_plot
        roc_metric_batches = []
        roc_metrics = results.roc_metrics.sort_values('ROC_AUC',
                                                      ascending=False)
        for idx in range(0, len(roc_metrics), num_curves_per_plot):
            batch = roc_metrics.iloc[idx:idx+num_curves_per_plot]
            roc_metric_batches.append(batch)
        roc_curves_file_name = None if dir is None else os.path.join(
            dir, "roc_curves_")
        num_batches = len(roc_metric_batches)
        for batch_no, batch in enumerate(roc_metric_batches):
            batch_file_name = None if roc_curves_file_name is None else \
                roc_curves_file_name + str(batch_no) + ".png"
            self._plot_roc_curves(batch, results.user_vep_name,
                                  results.roc_curve_coordinates,
                                  batch_no, num_batches,
                                  ves_color_palette, batch_file_name)
        roc_table_file_name = None if dir is None else os.path.join(
            dir, "roc_table" + ".html")
        self._display_roc_table(results, roc_table_file_name)

    def plot_mwu_results(self, results: VEAnalysisResult,
                         ves_color_palette: dict,
                         dir: str = None):
        num_bars_per_plot = self._mwu_config.num_bars_per_plot
        batches = []
        metrics = results.mwu_metrics.query('EXCEPTION.isna()')
        metrics = metrics.sort_values('NEG_LOG10_MWU_PVAL',
                                      ascending=False)
        for idx in range(0, len(metrics), num_bars_per_plot):
            batch = metrics.iloc[idx:idx+num_bars_per_plot]
            batches.append(batch)
        mwu_bar_file_name = None if dir is None else os.path.join(
            dir, "mwu_bar_")
        num_batches = len(batches)
        for batch_no, batch in enumerate(batches):
            batch_file_name = None if mwu_bar_file_name is None else \
                mwu_bar_file_name + str(batch_no) + ".png"
            self._plot_mwu_bar(batch, batch_no, num_batches, batch_file_name,
                               ves_color_palette)
        mwu_table_file_name = None if dir is None else os.path.join(
            dir, "mwu_table" + ".html")
        self._display_mwu_table(results, mwu_table_file_name)

    def plot_results(self, results: VEAnalysisResult,
                     metrics: str | list[str] = ["roc", "pr", "mwu"],
                     num_top_genes: int = 10,
                     dir: str = None):
        """
        Plot the results of an analysis either to the screen or to
        files.

        Parameters
        ----------
        results : VEAnalysisResult
            Analysis result object
        metrics : str or list[str]
            Specifies which metrics to plot. Can be a string
            indicating a single metric or a list of strings for
            multiple metrics. The metrics are: roc, pr, mwu.
        num_top_genes : int
            Number of top genes to plot. The top genes are the ones for
            which the most variants were observed.
        dir : str, optional
            Directory to place the plot files. The files will
            be placed in a subdirectory off of this directory
            whose name begins with ve_analysis_plots and suffixed
            by a unique timestamp. If not specified will plot
            to screen.
        """
        if type(metrics) is str:
            metrics = [metrics]
        if dir is not None:
            dir = unique_file_name(dir, "ve_analysis_plots_")
            create_folder(dir)
        ves_color_palette = create_feature_palette(
            results.general_metrics["SOURCE_NAME"])
        if "roc" in metrics and results.roc_metrics is not None:
            self.plot_roc_results(results, ves_color_palette, dir)
        if "pr" in metrics and results.pr_metrics is not None:
            self.plot_pr_results(results, ves_color_palette, dir)
        if "mwu" in metrics and results.mwu_metrics is not None:
            self.plot_mwu_results(results, ves_color_palette, dir)
        if results.gene_general_metrics is not None:
            self.plot_gene_results(results, metrics, num_top_genes, dir)

    def plot_gene_results(self, results: VEAnalysisResult,
                          metrics: list[str],
                          num_top_genes: int,
                          dir: str = None):
        """
        Plot gene-level results of an analysis.

        Parameters
        ----------
        results : VEAnalysisResult
            Analysis result object containing gene-level metrics
        metrics : list[str]
            List of metrics to plot (roc, pr, mwu)
        num_top_genes : int
            Number of top genes to plot
        dir : str, optional
            Directory to place the plot files
        """
        ves_color_palette = create_feature_palette(
            results.gene_general_metrics["SOURCE_NAME"].unique())
        gene_metric_sorter = GeneMetricSorter(
            results.gene_unique_variant_counts_df, num_top_genes)

        if "roc" in metrics and results.gene_roc_metrics is not None:
            self.plot_gene_roc_results(results, num_top_genes,
                                       gene_metric_sorter, 
                                       ves_color_palette, dir)

        if "pr" in metrics and results.gene_pr_metrics is not None:
            self.plot_gene_pr_results(results, gene_metric_sorter,
                                      ves_color_palette, dir)

        if "mwu" in metrics and results.gene_mwu_metrics is not None:
            self.plot_gene_mwu_results(results, gene_metric_sorter,
                                       ves_color_palette, dir)

    def plot_gene_roc_results(self, results: VEAnalysisResult,
                              num_top_genes: int,
                              gene_metric_sorter: GeneMetricSorter,
                              ves_color_palette: dict,
                              dir: str = None):

        gene_roc_metrics = results.gene_general_metrics.merge(
            results.gene_roc_metrics, how="inner", suffixes=(None, "_y"),
            on=['SCORE_SOURCE','GENE_SYMBOL'])[[
                'SOURCE_NAME', 'GENE_SYMBOL','NUM_VARIANTS', 'NUM_POSITIVE_LABELS',
                'NUM_NEGATIVE_LABELS', 'ROC_AUC']]
        gene_roc_metrics = gene_metric_sorter.sort_gene_metrics(
            gene_roc_metrics)

        heatmap_file_name = None if dir is None else os.path.join(
            dir, "gene_roc_heatmap.png")
        
        heatmap(gene_roc_metrics, 'SOURCE_NAME', 'GENE_SYMBOL', 'ROC_AUC',
                title='Gene-Level ROC', filename=heatmap_file_name)
                                          
        # Create gene ROC table
        table_file_name = None if dir is None else os.path.join(
            dir, "gene_roc_table.html")
        self._display_gene_roc_table(gene_roc_metrics, table_file_name)

    def plot_gene_pr_results(self, results: VEAnalysisResult,
                              ves_color_palette: dict,
                              dir: str = None):
        # Implement gene PR curve plotting
        pass

    def plot_gene_mwu_results(self, results: VEAnalysisResult,
                              ves_color_palette: dict,
                              dir: str = None):
        # Implement gene MWU bar plotting
        pass

    def _display_gene_roc_table(self, roc_results: pd.DataFrame,
                                file_name: str = None):
        style = roc_results.style.set_properties(
            subset=["NUM_VARIANTS",
                    "NUM_POSITIVE_LABELS", "NUM_NEGATIVE_LABELS",
                    "ROC_AUC"], **{"text-align": "right"})
        style = style.hide().relabel_index(
             ['VEP', "Gene", "Variant Total", "Positive Labels", "Negative Labels",
              "ROC AUC"],
             axis=1).set_caption("Gene-Level ROC")
        if file_name:
            style.to_html(file_name)
        else:
            display(style)


