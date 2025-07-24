import os

from .model import VEAnalysisCalibrationResult, VEAnalysisResult
from .file_util import (
    unique_file_name,
    create_folder
)


class VEAnalysisExporter:
    """Export results of an analysis to data files"""

    def export_results(self, results: VEAnalysisResult,
                       dir: str):
        """
        Export the results of an analysis to data files.

        Parameters
        ----------
        results : VEAnalysisResult
            Analysis result object with all relevant metrics
        dir : str
            Directory to place the data files. The files will
            be placed in a subdirectory off of this directory
            whose name begins with ve_analysis_data and suffixed
            by a unique timestamp.
        """
        dir = unique_file_name(dir, "ve_analysis_data_")
        create_folder(dir)
        results.general_metrics.to_csv(
            os.path.join(dir, "general_metrics.csv"), index=False)
        if results.roc_metrics is not None:
            results.roc_metrics.to_csv(
                os.path.join(dir, "roc_metrics.csv"), index=False)
            results.roc_curve_coordinates.to_csv(
                os.path.join(dir, "roc_curve_coords.csv"), index=False)
        if results.pr_metrics is not None:
            results.pr_metrics.to_csv(
                os.path.join(dir, "pr_metrics.csv"), index=False)
            results.pr_curve_coordinates.to_csv(
                os.path.join(dir, "pr_curve_coords.csv"), index=False)
        if results.mwu_metrics is not None:
            results.mwu_metrics.to_csv(
                os.path.join(dir, "mwu_metrics.csv"), index=False)
        if results.variants_included is not None:
            results.variants_included.to_csv(
                os.path.join(dir, "included_variants.csv"), index=False)
        if results.gene_general_metrics is not None:
            results.gene_general_metrics.to_csv(
                os.path.join(dir, "gene_general_metrics.csv"), index=False)
        if results.gene_roc_metrics is not None:
            results.gene_roc_metrics.to_csv(
                os.path.join(dir, "gene_roc_metrics.csv"), index=False)
        if results.gene_pr_metrics is not None:
            results.gene_pr_metrics.to_csv(
                os.path.join(dir, "gene_pr_metrics.csv"), index=False)
        if results.gene_mwu_metrics is not None:
            results.gene_mwu_metrics.to_csv(
                os.path.join(dir, "gene_mwu_metrics.csv"), index=False)
        if results.gene_roc_curve_coordinates is not None:
            results.gene_roc_curve_coordinates.to_csv(
                os.path.join(dir, "gene_roc_curve_coords.csv"), index=False)
        if results.gene_pr_curve_coordinates is not None:
            results.gene_pr_curve_coordinates.to_csv(
                os.path.join(dir, "gene_pr_curve_coords.csv"), index=False)

    def export_calibration_results(
            self, results: VEAnalysisCalibrationResult, dir: str):
        """
        Export the results of a calibration analysis to data files.

        Parameters
        ----------
        results : VEAnalysisCalibrationResult
            Calibration result object returned by calling
            VEAnalyzer.compute_calibration_metrics.
        dir : str
            Directory to place the data files. The files will
            be placed in a subdirectory off of this directory
            whose name begins with ve_calibration_data and suffixed
            by a unique timestamp.
        """
        dir = unique_file_name(dir, "ve_calibration_data_")
        create_folder(dir)
        results.pr_curve_coordinates_df.to_csv(
            os.path.join(dir, "pr_curve_coords.csv"), index=False)
        results.f1_curve_coordinates_df.to_csv(
            os.path.join(dir, "f1_curve_coords.csv"), index=False)
        results.score_pathogenic_fraction_df.to_csv(
            os.path.join(dir, "binned_scores.csv"),
            index=False)
        results.scores_and_labels_df.to_csv(
            os.path.join(dir, "included_variants.csv"),
            index=False)
