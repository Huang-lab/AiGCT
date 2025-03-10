import context  # noqa: F401

from aigct.etl.repo_loader import RepositoryLoader
from aigct.etl.container import VEETLContainer

LEGACY_DATA_FOLDERS = {
    "CANCER": "TGCA.V1/datas/cancer",
    "ADRD": "TGCA.V1/data_new/AD",
    "CHD": "TGCA.V1/data_new/CHD",
    "DDD": "TGCA.V1/data_new/DDD",
}
LEGACY_CANCER_VARIANT_FILES = [
    {"source": "HOTSPOT", "label": 1, "file": "hotspot.csv"},
    {"source": "MSK_PASSENGER", "label": 0, "file": "MSK_passenger.csv"},
    {"source": "TCGA_PASSENGER", "label": 0, "file": "TCGA_passenger.csv"},
]
LEGACY_ADRD_VARIANT_FILES = [
    {"source": "ADRD", "label": 1, "file": "ADcase.tsv"},
    {"source": "ADRD", "label": 0, "file": "ADcontrol.tsv"}
]
LEGACY_CHD_VARIANT_FILES = [
    {"source": "CHD", "label": 1, "file": "CHDcase.tsv"},
    {"source": "CHD", "label": 0, "file": "CHDcontrol.tsv"}
]
LEGACY_DDD_VARIANT_FILES = [
    {"source": "DDD", "label": 1, "file": "DDDcase.tsv"},
    {"source": "DDD", "label": 0, "file": "DDDcontrol.tsv"}
]

LEGACY_VARIANT_FILES = {
    "CANCER": LEGACY_CANCER_VARIANT_FILES,
    "ADRD": LEGACY_ADRD_VARIANT_FILES,
    "CHD": LEGACY_CHD_VARIANT_FILES,
    "DDD": LEGACY_DDD_VARIANT_FILES
}


def migrate_task_files(loader: RepositoryLoader, task: str):
    for file in LEGACY_VARIANT_FILES[task]:
        loader.load_variant_file("hg38", task, file["file"],
                                 LEGACY_DATA_FOLDERS[task],
                                 file["source"], file["label"],
                                 "hg19", "hg18")


container = VEETLContainer()
loader = container.loader
"""
loader.init_variant_task()
loader.init_variant_effect_source()
migrate_task_files(loader, "CANCER")
migrate_task_files(loader, "ADRD")
migrate_task_files(loader, "CHD")
migrate_task_files(loader, "DDD")
"""
task = "CANCER"
loader.load_variant_file("hg38", task, "hotspot.csv",
                         LEGACY_DATA_FOLDERS[task],
                         "hotspot", 1,
                         "hg19", "hg18")
loader.load_variant_file("hg38", task, "MSK_passenger.csv",
                         LEGACY_DATA_FOLDERS[task],
                         "MSK_PASSENGER", 0,
                         "hg19", "hg18")
