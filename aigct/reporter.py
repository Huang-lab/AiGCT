import pandas as pd
from .model import (
    VEAnalysisResult
)
from .repository import (
     VariantEffectScoreRepository,
     VariantEffectLabelRepository,
     VariantEffectSourceRepository,
     VARIANT_PK_COLUMNS
)
from .date_util import now_str_basic_format, now_str_compact
import os
import sys
from .file_util import new_line

VARIANT_EFFECT_SCORE_COLS = ["SCORE_SOURCE"] +\
    VARIANT_PK_COLUMNS + ["RANK_SCORE"]


class VEAnalysisReporter:
    """Report analysis results"""

    def _write_metric_dataframe(self, out, metric_df: pd.DataFrame):
        out.write(metric_df.to_string(index=False))
        new_line(out)

    def _write_summary(self, out, metrics: VEAnalysisResult):
        new_line(out)
        out.write("Summary metrics for Variant Effect Prediction Benchmark: " +
                  now_str_basic_format())
        new_line(out, 2)
        out.write("Total number of user supplied variants: " +
                  str(metrics.num_user_variants))
        new_line(out, 2)
        out.write("Total number of variants across all VEPs in analysis: " +
                  str(metrics.num_variants_included))
        new_line(out, 2)
        self._write_metric_dataframe(out, metrics.general_metrics)
        new_line(out, 2)
        if metrics.roc_metrics is not None:
            out.write("ROC Metrics")
            new_line(out, 2)
            self._write_metric_dataframe(
                out, metrics.roc_metrics.sort_values(by='ROC_AUC',
                                                     ascending=False))
            new_line(out, 2)
        if metrics.pr_metrics is not None:
            out.write("Precision/Recall Metrics")
            new_line(out, 2)
            self._write_metric_dataframe(
                out, metrics.pr_metrics.sort_values(by="PR_AUC",
                                                    ascending=False))
            new_line(out, 2)
        if metrics.mwu_metrics is not None:
            out.write("Mann-Whitney U -log10(P value)")
            new_line(out, 2)
            self._write_metric_dataframe(
                out, metrics.mwu_metrics.sort_values(by="NEG_LOG10_MWU_PVAL",
                                                     ascending=False))
            new_line(out, 2)

    def write_summary(self, metrics: VEAnalysisResult,
                      dir: str = None):
        """
        Generate a report summarizing the results of an analysis. It will
        be written either to the screen or to a file.

        Parameters
        ----------
        results : VEAnalysisResult
            Object containing the results of an analysis.
        dir : str, optional
            Directory to place the report file. The file name will
            begin with variant_bm_summary and suffixed
            by a unique timestamp. If not specified will print to the
            screen.
        """
        if dir is not None:
            outfile = os.path.join(dir, now_str_compact("variant_bm_summary")
                                   + ".txt")
            with open(outfile, "w") as out:
                self._write_summary(out, metrics)
        else:
            self._write_summary(sys.stdout, metrics)
        
        
        


