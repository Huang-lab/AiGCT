"""Driver for the benchmark curation pipeline.

Reads datasets.yaml and runs the requested stage. Each stage is idempotent and
writes into the output directories configured in config.yaml, so stages can be
re-run individually.

    python -m curation.run_curation annotate            # all cohort datasets
    python -m curation.run_curation annotate ASD_case1  # just one
    python -m curation.run_curation clinvar             # scan dbNSFP for ClinVar
    python -m curation.run_curation balance             # gene-balance ClinVar
    python -m curation.run_curation overlap             # positive/negative split
    python -m curation.run_curation exclude-clinvar     # Figure 4B/D/F inputs
    python -m curation.run_curation all

Note that `annotate` and `clinvar` both read the full dbNSFP release and take
hours; `clinvar` scans all 25 chromosome files start to finish.
"""
import argparse
import os

import yaml

from . import balance, clinvar, dbnsfp, overlap
from .config import REPO_ROOT, annotation_path, output_path, repo_path

DATASETS_FILE = os.path.join(REPO_ROOT, "datasets.yaml")


def load_datasets():
    with open(DATASETS_FILE, encoding="utf-8") as f:
        return yaml.safe_load(f)


def stage_annotate(spec, only=None):
    """Annotate cohort coordinate lists against dbNSFP."""
    for name, entry in spec["datasets"].items():
        if "annotation" not in entry:
            continue  # ClinVar datasets are not built from a coordinate list
        if only and name not in only:
            continue
        print(f"\n=== {name} ({entry['assembly']}) ===")
        dbnsfp.annotate_dataset(
            annotation_file=annotation_path(entry["annotation"]),
            out_path=output_path("processed", entry["output"]),
            assembly=entry["assembly"],
            dataset_name=name,
        )


def stage_clinvar(spec):
    """Scan the dbNSFP release for ClinVar-classified variants."""
    entry = spec["datasets"]["clinvar"]
    clinvar.build(out_path=output_path("processed", entry["output"]))


def stage_balance(spec):
    """Sample an equal number of pathogenic and benign variants per gene."""
    import pandas as pd

    entry = spec["datasets"]["balanced_clinvar"]
    df = pd.read_csv(output_path("processed", entry["input"]))
    balance.balance_by_gene(
        df,
        output=output_path("processed", entry["output"]),
        gene_stats_path=output_path("processed", "gene_variant_counts.csv"),
        excluded_genes_path=output_path("processed", "excluded_genes.csv"),
    )


def stage_overlap(spec):
    """Drop variants shared between the positive and negative sets of a task."""
    for job in spec["overlap_removal"]:
        print(f"\n=== overlap removal: {job['task']} ===")
        overlap.remove_overlapping(job["positives"], job["negatives"])


def stage_exclude_clinvar(spec):
    """Build the ClinVar-excluded developmental disorder benchmarks."""
    cfg = spec["clinvar_exclusion"]
    for task, entry in cfg["datasets"].items():
        print(f"\n=== ClinVar exclusion: {task} ===")
        overlap.exclude_clinvar_variants(
            dataset_files=entry["inputs"],
            clinvar_files=cfg["clinvar_set"],
            output_name=entry["output"],
        )


STAGES = {
    "annotate": lambda spec, only: stage_annotate(spec, only),
    "clinvar": lambda spec, only: stage_clinvar(spec),
    "balance": lambda spec, only: stage_balance(spec),
    "overlap": lambda spec, only: stage_overlap(spec),
    "exclude-clinvar": lambda spec, only: stage_exclude_clinvar(spec),
}

ALL_ORDER = ["annotate", "clinvar", "balance", "overlap", "exclude-clinvar"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=list(STAGES) + ["all"])
    parser.add_argument(
        "datasets",
        nargs="*",
        help="restrict the annotate stage to these dataset names",
    )
    args = parser.parse_args()

    spec = load_datasets()
    stages = ALL_ORDER if args.stage == "all" else [args.stage]
    for stage in stages:
        print(f"\n########## {stage} ##########")
        STAGES[stage](spec, args.datasets or None)


if __name__ == "__main__":
    main()
