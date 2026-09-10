"""Configuration loading and path resolution.

All modules in this package read their paths from config.yaml rather than
hard-coding them, so the pipeline can be run from any checkout location.
Set AIGCT_CURATION_CONFIG to point at an alternative config file.
"""
import os
from functools import lru_cache

import yaml

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_CONFIG = os.path.join(REPO_ROOT, "config.yaml")


@lru_cache(maxsize=None)
def load_config(path=None):
    path = path or os.environ.get("AIGCT_CURATION_CONFIG", DEFAULT_CONFIG)
    with open(path, encoding="utf-8") as f:
        return yaml.safe_load(f)


def repo_path(*parts):
    """Resolve a path relative to the repository root."""
    return os.path.join(REPO_ROOT, *parts)


def annotation_path(name):
    cfg = load_config()
    return repo_path(cfg["repo"]["annotation_dir"], name)


def output_path(kind, name=None):
    """Resolve a path under one of the configured output directories.

    `kind` is one of the keys under `output:` in config.yaml, without the
    trailing "_dir" (e.g. "processed", "duplicates", "logs").
    """
    cfg = load_config()
    directory = repo_path(cfg["output"][kind + "_dir"])
    os.makedirs(directory, exist_ok=True)
    return directory if name is None else os.path.join(directory, name)


def dbnsfp_chromosome_file(chrom):
    """Path to the dbNSFP release file for one chromosome."""
    cfg = load_config()["external"]
    return os.path.join(
        cfg["dbnsfp_dir"],
        cfg["dbnsfp_file_template"].format(chrom=chrom),
    )
