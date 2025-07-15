import matplotlib.pyplot as plt
# from IPython.display import HTML
from IPython.display import display
# import seaborn as sns
# import plotly.express as px
from matplotlib.ticker import MaxNLocator
import pandas as pd
import os
from .model import (
    VEAnalysisCalibrationResult,
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
                         num_top_labelled_veps: int,
                         ves_color_palette: dict,
                         file_name: str = None):
        title = self._roc_pr_config.roc_title + (
            f" ({batch_no+1} of {num_batches})" if num_batches > 1 else "")
        plt.figure(figsize=(10, 8))
        plot_count = 0
        for i, ve_auc in aucs.iterrows():
            if pd.isna(ve_auc['ROC_AUC']):
                continue
            ve_curve_coords = curve_coords[curve_coords['SCORE_SOURCE'] ==
                                           ve_auc['SCORE_SOURCE']].sort_values(
                'THRESHOLD', ascending=False)
            if num_top_labelled_veps is None or plot_count < num_top_labelled_veps:
                label = ve_auc['SOURCE_NAME'] + \
                    ' (AUC=' + str(round(ve_auc['ROC_AUC'], 4)) + ')'
            else:
                label = None
            plt.plot(
                ve_curve_coords['FALSE_POSITIVE_RATE'],
                ve_curve_coords['TRUE_POSITIVE_RATE'],
                label=label,
                color=ves_color_palette[ve_auc['SOURCE_NAME']],
                lw=self._config.line_width,
                linestyle=self._config.line_style)
            plot_count += 1

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
                        num_top_labelled_veps: int,
                        ves_color_palette: dict, file_name: str = None):
        title = self._roc_pr_config.pr_title + \
            (f" ({batch_no+1} of {num_batches})" if num_batches > 1 else "")
        plt.figure(figsize=(10, 8))
        plot_count = 0
        for i, ve_auc in aucs.iterrows():
            ve_curve_coords = curve_coords[curve_coords['SCORE_SOURCE'] ==
                                           ve_auc['SCORE_SOURCE']].sort_values(
                'THRESHOLD', ascending=False)
            if num_top_labelled_veps is None or plot_count < num_top_labelled_veps:
                label = ve_auc['SOURCE_NAME'] + \
                    ' (AUC=' + str(round(ve_auc['PR_AUC'], 4)) + ')'
            else:
                label = None
            plt.plot(
                ve_curve_coords['RECALL'][1:],
                ve_curve_coords['PRECISION'][1:],
                label=label,
                color=ves_color_palette[ve_auc['SOURCE_NAME']],
                lw=self._config.line_width,
                linestyle=self._config.line_style)
            plot_count += 1

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
                'NUM_NEGATIVE_LABELS', 'NEG_LOG10_MWU_PVAL']]
        table_output = table_output.sort_values('NEG_LOG10_MWU_PVAL',
                                                ascending=False)
        style = table_output.style.set_properties(
            subset=["NUM_VARIANTS",
                    "NUM_POSITIVE_LABELS", "NUM_NEGATIVE_LABELS",
                    "NEG_LOG10_MWU_PVAL"], **{"text-align": "right"})
        style = style.hide().relabel_index(
             ['VEP', "Variant Total", "Positive Labels", "Negative Labels",
              "MWU -log10(pval)"],
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
                 y_tick_label_size=config.y_tick_label_size,
                 x_tick_label_size=config.x_tick_label_size,
                 xtick_rotation=config.xtick_rotation,
                 xtick_rotation_mode=config.xtick_rotation_mode,
                 file_dpi=self._config.file_dpi,
                 bbox_inches=self._config.bbox_inches)

    def plot_pr_results(self, results: VEAnalysisResult,
                        num_top_labelled_veps: int,
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
                                 num_top_labelled_veps,
                                 ves_color_palette, batch_file_name)
        pr_table_file_name = None if dir is None else os.path.join(
            dir, "pr_table" + ".html")
        self._display_pr_table(results, pr_table_file_name)

    def plot_roc_results(self, results: VEAnalysisResult,
                         num_top_labelled_veps: int,
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
                                  num_top_labelled_veps,
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
                     num_top_labelled_veps: int = None,
                     num_top_genes: int = None,
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
        num_top_labelled_veps : int
            If not None, only this many of the top performing veps will be
            labelled in the auc plot legends. This is useful when there are
            many veps and the legend becomes too cluttered.
        num_top_genes : int
            If compute_gene_metrics was set to True in call to compute_metrics,
            then only include this many top genes in the plot. The top gene
            are the ones for which the most variants were observed.
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
            self.plot_roc_results(results, num_top_labelled_veps,
                                  ves_color_palette, dir)
        if "pr" in metrics and results.pr_metrics is not None:
            self.plot_pr_results(results, num_top_labelled_veps,
                                 ves_color_palette, dir)
        if "mwu" in metrics and results.mwu_metrics is not None:
            self.plot_mwu_results(results, ves_color_palette, dir)
        if results.gene_general_metrics is not None:
            self.plot_gene_results(results, metrics, num_top_genes, dir)

    def plot_gene_results(self, results: VEAnalysisResult,
                          metrics: list[str],
                          num_top_genes: int = None,
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
            Number of top genes to plot based on the number of variants
            in each gene included in the analysis.
        dir : str, optional
            Directory to place the plot files
        """
        ves_color_palette = create_feature_palette(
            results.gene_general_metrics["SOURCE_NAME"].unique())
        gene_metric_sorter = GeneMetricSorter(
            results.gene_unique_variant_counts_df, num_top_genes)

        if "roc" in metrics and results.gene_roc_metrics is not None:
            self.plot_gene_level_results(results.gene_general_metrics,
                                         results.gene_roc_metrics,
                                         "ROC_AUC",
                                         "ROC AUC",
                                         "Gene-Level ROC",
                                         gene_metric_sorter,
                                         ves_color_palette,
                                         dir,
                                         "gene_roc_heatmap.png",
                                         "gene_roc_table.html")

        if "pr" in metrics and results.gene_pr_metrics is not None:
            self.plot_gene_level_results(results.gene_general_metrics,
                                         results.gene_pr_metrics,
                                         "PR_AUC",
                                         "PR AUC",
                                         "Gene-Level Precision/Recall",
                                         gene_metric_sorter,
                                         ves_color_palette,
                                         dir,
                                         "gene_pr_heatmap.png",
                                         "gene_pr_table.html")

        if "mwu" in metrics and results.gene_mwu_metrics is not None:
            self.plot_gene_level_results(results.gene_general_metrics,
                                         results.gene_mwu_metrics,
                                         "NEG_LOG10_MWU_PVAL",
                                         "MWU -log10(pval)",
                                         "Gene-Level Mann-Whitney U",
                                         gene_metric_sorter,
                                         ves_color_palette,
                                         dir,
                                         "gene_mwu_heatmap.png",
                                         "gene_mwu_table.html")

    def plot_gene_level_results(self, gene_general_metrics: pd.DataFrame,
                                gene_metrics: pd.DataFrame,
                                metric_column: str,
                                metric_display_name: str,
                                title: str,
                                gene_metric_sorter: GeneMetricSorter,
                                ves_color_palette: dict,
                                dir: str = None,
                                figure_file_name: str = None,
                                table_file_name: str = None):
        gene_metrics = gene_general_metrics.merge(
            gene_metrics, how="inner", suffixes=(None, "_y"),
            on=['SCORE_SOURCE', 'GENE_SYMBOL'])[[
                'SOURCE_NAME', 'GENE_SYMBOL', 'NUM_VARIANTS',
                'NUM_POSITIVE_LABELS',
                'NUM_NEGATIVE_LABELS', metric_column]]
        gene_metrics = gene_metric_sorter.sort_gene_metrics(
            gene_metrics)

        heatmap_file_name = None if dir is None else os.path.join(
            dir, figure_file_name)

        gene_metrics_heatmap = gene_metrics.query(f'{metric_column}.notna()')
        heatmap(gene_metrics_heatmap, 'SOURCE_NAME', 'GENE_SYMBOL',
                metric_column,
                title=title, x_label='Gene', y_label='VEP',
                cbar_kws={'label': metric_display_name,
                          'shrink': 0.8},
                filename=heatmap_file_name)

        # Create gene ROC table
        table_file_name = None if dir is None else os.path.join(
            dir, table_file_name)
        self._display_gene_metric_table(gene_metrics, metric_column,
                                        metric_display_name, title,
                                        table_file_name)

    def _display_gene_metric_table(self, metrics_df: pd.DataFrame,
                                   metric_column: str,
                                   metric_display_name: str,
                                   title: str,
                                   file_name: str = None):
        style = metrics_df.style.set_properties(
            subset=["NUM_VARIANTS",
                    "NUM_POSITIVE_LABELS", "NUM_NEGATIVE_LABELS",
                    metric_column], **{"text-align": "right"})
        style = style.hide().relabel_index(
             ['VEP', "Gene", "Variant Total", "Positive Labels", "Negative Labels",
              metric_display_name],
             axis=1).set_caption(title)
        if file_name:
            style.to_html(file_name)
        else:
            display(style)

    def _plot_score_vs_pathogenic_fraction(self, results, annotate: bool = False,
                                           dir: str = None):
        """
        Plot the fraction of pathogenic variants versus mean score as a line
        plot.

        Parameters
        ----------
        results : VEAnalysisCalibrationResult
            Calibration result object containing score bins and pathogenic
            fractions
        dir : str, optional
            Directory to save the plot file
        """
        title = "Pathogenic fraction by score"
        sorted_df = results.score_pathogenic_fraction_df.sort_values(
            "MEAN_SCORE")
        pathogenic_fraction = (
            sorted_df["NUM_POSITIVE_LABELS"] / sorted_df["NUM_VARIANTS"]
            ).round(1)
        mean_score = sorted_df["MEAN_SCORE"].round(2)

        # Create the plot
        plt.figure(figsize=(10, 6))
        
        # Plot pathogenic fraction vs mean score as a line plot
        plt.plot(mean_score, pathogenic_fraction, 
                 marker='o', linestyle='-', linewidth=2, 
                 markersize=6, color='blue')

        if annotate:
            for row in sorted_df.itertuples():
                plt.annotate(row.SCORE_RANGE,
                             (row.MEAN_SCORE,
                              row.NUM_POSITIVE_LABELS/row.NUM_VARIANTS),
                             xytext=(5, 5),  # Offset from point
                             textcoords='offset points',
                             fontsize=8,
                             ha='left',
                             va='bottom',
                             bbox=dict(boxstyle='round,pad=0.3',
                                       facecolor='yellow',
                                       alpha=0.7),
                             rotation=45)  # Rotate text to avoid overlap
        
        # Add diagonal reference line (y=x)
        x_min, x_max = plt.xlim()
        y_min, y_max = plt.ylim()
        reference_line_x = [x_min, x_max]
        reference_line_y = [y_min, y_max]
        plt.plot(reference_line_x, reference_line_y, linestyle='--',
                 color='gray', alpha=0.7)

        # Set labels and title
        plt.xlabel('Score', fontsize=14)
        plt.ylabel('Pathogenic Fraction', fontsize=14)
        plt.title(title, fontsize=16)
        
        # Set y-axis limits
        plt.ylim(0, 1.05)
        
        # Add grid
        plt.grid(True, linestyle='--', alpha=0.7)
        
        # Save or display the plot
        if dir is None:
            plt.show()
        else:
            file_path = os.path.join(dir, "pathogenic_fraction_by_score.png")
            plt.savefig(file_path, dpi=self._config.file_dpi,
                        bbox_inches=self._config.bbox_inches)
            plt.savefig(file_path.replace(".png", ".svg"),
                        format='svg',
                        bbox_inches=self._config.bbox_inches)
        
        plt.close()

    def _plot_precision_recall_vs_thresholds1(
            self, results: VEAnalysisCalibrationResult, dir: str = None):
        """
        Plot precision, recall, and F1 score versus threshold values.
        
        Parameters
        ----------
        results : VEAnalysisCalibrationResult
            Calibration result object containing precision-recall curve data
        dir : str, optional
            Directory to save the plot file
        """
        # Get precision-recall curve data
        pr_curve_coords = results.negative_pr_curve_coordinates
        
        # Sort by threshold
        pr_curve_coords = pr_curve_coords.sort_values('THRESHOLD')
        
        # Calculate F1 score
        precision = pr_curve_coords['PRECISION']
        recall = pr_curve_coords['RECALL']
        f1_scores = 2 * (precision * recall) / (precision + recall + 1e-8)  # avoid division by zero
        
        # Create the plot
        plt.figure(figsize=(10, 6))
        
        # Plot precision, recall, and F1 score vs threshold
        plt.plot(pr_curve_coords['THRESHOLD'], precision, 
                 label='Precision', color='blue', linewidth=2)
        plt.plot(pr_curve_coords['THRESHOLD'], recall, 
                 label='Recall', color='green', linewidth=2)
        plt.plot(pr_curve_coords['THRESHOLD'], f1_scores, 
                 label='F1 Score', color='red', linewidth=2)
        
        # Set labels and title
        plt.xlabel('Decision Threshold', fontsize=14)
        plt.ylabel('Score', fontsize=14)
        plt.title('Performance Metrics vs. Decision Threshold', fontsize=16)
        
        # Add legend
        plt.legend(fontsize=12)
        
        # Add grid
        plt.grid(True, linestyle='--', alpha=0.7)
        
        # Set y-axis limits
        plt.ylim(0, 1.05)
        
        # Save or display the plot
        if dir is None:
            plt.show()
        else:
            file_path = os.path.join(dir, "precision_recall_vs_threshold.png")
            plt.savefig(file_path, dpi=self._config.file_dpi, 
                        bbox_inches=self._config.bbox_inches)
            plt.savefig(file_path.replace(".png", ".svg"), 
                        format='svg', 
                        bbox_inches=self._config.bbox_inches)
        
        plt.close()

    def _plot_precision_recall_vs_threshold(
            self, results: VEAnalysisCalibrationResult,
            threshold_boundary: float = 0.5,
            precision_cutoff: float = 0.9,
            dir: str = None):
        """
        Plot precision, recall, and F1 score versus threshold values.
        
        Parameters
        ----------
        results : VEAnalysisCalibrationResult
            Calibration result object containing precision-recall curve data
        dir : str, optional
            Directory to save the plot file
        """
        # Get precision-recall curve data
        pr_curve_coords = results.negative_pr_curve_coordinates.query(
                "THRESHOLD <= @threshold_boundary")
        
        # Sort by threshold
        pr_curve_coords = pr_curve_coords.sort_values('THRESHOLD')
        
        # Calculate F1 score
        precision = pr_curve_coords['PRECISION']
        recall = pr_curve_coords['RECALL']
        
        # Create the plot
        plt.figure(figsize=(10, 6))
        
        # Plot precision, recall, and F1 score vs threshold
        plt.plot(pr_curve_coords['THRESHOLD'], precision, 
                 label='Benign Precision', color='blue', linewidth=2)
        plt.plot(pr_curve_coords['THRESHOLD'], recall, 
                 label='Benign Recall', color='blue', linewidth=2,
                 linestyle='dashdot')
        if len(pr_curve_coords) > 1 and precision.max() >= precision_cutoff \
                and precision.min() <= precision_cutoff:
            plt.axvline(x=pr_curve_coords[precision >= precision_cutoff]['THRESHOLD'].min(),
                        color='gray', linestyle='--', linewidth=1, alpha=0.7)
        
        pr_curve_coords = results.positive_pr_curve_coordinates.query(
                "THRESHOLD >= @threshold_boundary")
        
        # Sort by threshold
        pr_curve_coords = pr_curve_coords.sort_values('THRESHOLD')
        
        # Calculate F1 score
        precision = pr_curve_coords['PRECISION']
        recall = pr_curve_coords['RECALL']
        
        # Plot precision, recall, and F1 score vs threshold
        plt.plot(pr_curve_coords['THRESHOLD'], precision, 
                 label='Pathogenic Precision', color='red', linewidth=2)
        plt.plot(pr_curve_coords['THRESHOLD'], recall, 
                 label='Pathogenic Recall', color='red', linewidth=2,
                 linestyle='dashdot')
        
        if len(pr_curve_coords) > 1 and precision.max() >= precision_cutoff \
                and precision.min() <= precision_cutoff:
            plt.axvline(x=pr_curve_coords[precision >= precision_cutoff]['THRESHOLD'].min(),
                        color='gray', linestyle='--', linewidth=1, alpha=0.7)

        plt.axhline(y=precision_cutoff,
                    color='gray', linestyle='--', linewidth=1, alpha=0.7)
        
        # Set labels and title
        plt.xlabel('Decision Threshold', fontsize=14)
        plt.ylabel('Precision or Recall', fontsize=14)
        plt.title('Performance Metrics vs. Decision Threshold', fontsize=16)
        
        # Add legend
        plt.legend(fontsize=12)
        
        # Set y-axis limits
        plt.ylim(0, 1.05)
        
        # Save or display the plot
        if dir is None:
            plt.show()
        else:
            file_path = os.path.join(dir, "precision_recall_vs_threshold.png")
            plt.savefig(file_path, dpi=self._config.file_dpi, 
                        bbox_inches=self._config.bbox_inches)
            plt.savefig(file_path.replace(".png", ".svg"), 
                        format='svg', 
                        bbox_inches=self._config.bbox_inches)
        
        plt.close()

    def _plot_score_vs_variant_counts(self,
                                      results: VEAnalysisCalibrationResult,
                                      bins: int,
                                      dir: str = None):
        """
        Plot true histograms showing the distribution of RANK_SCORE values
        for positive (BINARY_LABEL=1) and negative (BINARY_LABEL=0) variants.
        
        Parameters
        ----------
        results : VEAnalysisCalibrationResult
            Calibration result object containing variant scores and labels
        dir : str, optional
            Directory to save the plot file
        """
        # Get the raw data with individual variant scores
        df = results.scores_and_labels_df
        
        # Separate positive and negative variants
        positive_variants = df[df['BINARY_LABEL'] == 1]['RANK_SCORE']
        negative_variants = df[df['BINARY_LABEL'] == 0]['RANK_SCORE']
        
        # Create figure
        plt.figure(figsize=(10, 6))
        ax = plt.gca()
        
        # Define common histogram parameters
        alpha = 0.6
        
        # Plot histograms on the same axis with counts (not density)
        plt.hist(negative_variants, bins=bins, alpha=alpha, color='blue',
                 label=f'Benign variants (n={len(negative_variants)})')
        plt.hist(positive_variants, bins=bins, alpha=alpha, color='red',
                 label=f'Pathogenic variants (n={len(positive_variants)})')

        # Force y-axis ticks to be integers
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        
        # Add formatting
        plt.xlabel('Variant Score', fontsize=14)
        plt.ylabel('Variant Count', fontsize=14)
        plt.title('Distribution of Variant Scores by Pathogenicity',
                  fontsize=16)
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.legend(fontsize=12, loc='upper left')
        
        # Save or display the plot
        if dir is None:
            plt.show()
        else:
            file_path = os.path.join(dir, "score_variant_distribution.png")
            plt.savefig(file_path, dpi=self._config.file_dpi,
                        bbox_inches=self._config.bbox_inches)
            plt.savefig(file_path.replace(".png", ".svg"),
                        format='svg', 
                        bbox_inches=self._config.bbox_inches)
        
        plt.close()

    def plot_calibration_curves(self, results: VEAnalysisCalibrationResult,
                                pathogenic_benign_variant_count_bins = 30,
                                pr_threshold_boundary = 0.5,
                                precision_cutoff = 0.9,
                                dir: str = None):
        """
        Parameters
        ----------
        results : VEAnalysisResult
            Analysis result object
        variant_effect_source: str
            Specifies either the user vep name if the user supplied user
            vep scores in the call to VEAnalyzer.compute_metrics or the
            name of a system vep to plot.
        dir : str, optional
            Directory to place the plot files. The files will
            be placed in a subdirectory off of this directory
            whose name begins with ve_analysis_plots and suffixed
            by a unique timestamp. If not specified will plot
            to screen.
        """

        if dir is not None:
            dir = unique_file_name(dir, "ve_calibration_plots_")
            create_folder(dir)
        self._plot_score_vs_pathogenic_fraction(results,
                                                dir=dir)
        self._plot_score_vs_variant_counts(
            results, pathogenic_benign_variant_count_bins, dir
            )
        self._plot_metrics_vs_threshold(results, dir=dir)

    def _plot_metrics_vs_threshold(
            self, results: VEAnalysisCalibrationResult,
            dir: str = None):
        """
        Plot precision, recall, F1 score, and ROC metrics versus threshold values.
        
        Parameters
        ----------
        results : VEAnalysisCalibrationResult
            Calibration result object containing precision-recall and ROC curve data
        dir : str, optional
            Directory to save the plot file
        """
        # Create the plot
        plt.figure(figsize=(12, 8))

        # Plot Precision vs Threshold
        pr_coords = results.pr_curve_coordinates.sort_values('THRESHOLD')[:-1]
        plt.plot(pr_coords['THRESHOLD'], pr_coords['PRECISION'], 
                 label='Precision', color='blue', linewidth=2)
        plt.plot(pr_coords['THRESHOLD'], pr_coords['RECALL'], 
                 label='Recall', color='green', linewidth=2)

        roc_coords = results.roc_curve_coordinates.sort_values(
            'THRESHOLD')[:-1]
        plt.plot(roc_coords['THRESHOLD'], roc_coords['TRUE_POSITIVE_RATE'], 
                 label='True Positive Rate (ROC)', color='orange', linewidth=2)
        plt.plot(roc_coords['THRESHOLD'], roc_coords['FALSE_POSITIVE_RATE'], 
                 label='False Positive Rate (ROC)', color='red', linewidth=2)

        # Plot F1 Score vs Threshold
        f1_coords = results.f1_curve_coordinates.sort_values('THRESHOLD')[:-1]
        plt.plot(f1_coords['THRESHOLD'], f1_coords['F1_SCORE'], 
                 label='F1 Score', color='purple', linewidth=2)

        # Set labels and title
        plt.xlabel('Decision Threshold', fontsize=14)
        plt.ylabel('Metric Value', fontsize=14)
        plt.title('Performance Metrics vs. Decision Threshold', fontsize=16)

        # Add legend
        plt.legend(fontsize=12, loc='best')

        # Set y-axis limits
        plt.ylim(0, 1.05)

        # Save or display the plot
        if dir is None:
            plt.show()
        else:
            file_path = os.path.join(dir, "metrics_vs_threshold.png")
            plt.savefig(file_path, dpi=self._config.file_dpi, 
                        bbox_inches=self._config.bbox_inches)
            plt.savefig(file_path.replace(".png", ".svg"), 
                        format='svg', 
                        bbox_inches=self._config.bbox_inches)
        
        plt.close()
