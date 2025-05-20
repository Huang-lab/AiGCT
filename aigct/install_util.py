# import context
import requests
import re

import yaml
import random
import tarfile
import importlib.resources as pkg_resources
import io
import os
from .file_util import create_folder, unique_file_name
from .container import VEBenchmarkContainer
import aigct.config


def sample_user_scores(ve_bm_container):
    user_scores_df = ve_bm_container._score_repo.get(
        "CANCER", "REVEL")
    random_idxs = random.sample(list(range(len(user_scores_df))), 2000)
    user_variants = user_scores_df.iloc[random_idxs]
    user_variants['RANK_SCORE'] = user_variants['RANK_SCORE'].apply(
        lambda scor: scor + (random.uniform(0.05, 0.15) *
                             random.sample([1, -1], 1)[0])
        if scor < 0.84 and scor > 0.16 else scor)
    return user_variants


def get_sample_config() -> dict:
    with pkg_resources.open_text(aigct.config, 'aigct.yaml.sample') as sample:
    # with open("./config/aigct.yaml.sample") as sample:
        return yaml.safe_load(sample)


def init_config_file(config_dir: str = ".", data_dir: str = ".",
                     output_dir: str = "./analysis_output",
                     log_dir: str = "./log"):
    # create_folder(config_dir)
    config_file = os.path.join(config_dir, "aigct.yaml")
    config = get_sample_config()
    config["repository"]["root_dir"] = data_dir
    config["output_dir"] = output_dir
    config["log"]["dir"] = log_dir
    with (open(config_file, "w") as
            conf_file):
        yaml.dump(config, conf_file)


def init_db(conf_dir: str = "./config"):
    container = VEBenchmarkContainer(os.path.join(conf_dir, "aigct.yaml"))
    config = container.config
    url = config.repository.source_url
    version = config.repository.version
    dir = config.repository.root_dir
    response = requests.get(url, stream=True)
    try:
        header = response.headers["content-disposition"]
        p = re.compile(".+filename=((.+).tar.gz)")
        tar_file = p.match(header).group(1)
    except Exception:
        tar_file = "repo_" + version.replace(".", "_") + ".tar.gz"
    tar_file = os.path.join(dir, tar_file)
    with open(tar_file, "wb") as file:
        for chunk in response.iter_content(chunk_size=1024):
            if chunk:
                file.write(chunk)
    with tarfile.open(tar_file) as archive:
        archive.extractall(dir)


def check_install(conf_dir: str = "./config"):
    container = VEBenchmarkContainer(os.path.join(conf_dir, "aigct.yaml"))
    outdir = container.config.output_dir
    user_test_vep_scores = sample_user_scores(container)
    analyzer = container.analyzer
    metrics = analyzer.compute_metrics(
                "CANCER", user_test_vep_scores, vep_min_overlap_percent=50,
                variant_vep_retention_percent=1, list_variants=True)

    container.reporter.write_summary(metrics)
    container.plotter.plot_results(metrics)
    container.exporter.export_results(metrics, outdir)








