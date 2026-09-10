"""Filesystem locations used by the analysis notebooks.

Every location is read from an environment variable so the notebooks run on a
machine other than the one they were written on.  The defaults are relative to
this repository; point the variables at your own AIGCT installation instead.

    AIGCT_HOME         AIGCT installation that has already downloaded the
                       benchmark database (the directory holding config/ and db/)
    AIGCT_CONFIG       aigct.yaml of that installation
    AIGCT_DB_DIR       benchmark database tables (db/data/<TASK>/...)
    AIGCT_CLINVAR_CSV  ClinVar table from `python -m curation.run_curation clinvar`
    AIGCT_SCORES_DIR   MAVEN / EVE score tables
    AIGCT_OUTPUT_DIR   where the notebooks write their own output
"""
import os

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

AIGCT_HOME = os.environ.get("AIGCT_HOME", REPO_ROOT)

CONFIG = os.environ.get("AIGCT_CONFIG", os.path.join(AIGCT_HOME, "config", "aigct.yaml"))

DB_DIR = os.environ.get("AIGCT_DB_DIR", os.path.join(AIGCT_HOME, "db", "data"))

CLINVAR_CSV = os.environ.get(
    "AIGCT_CLINVAR_CSV",
    os.path.join(AIGCT_HOME, "output", "processed", "clinvar_withoutX.csv"),
)

SCORES_DIR = os.environ.get("AIGCT_SCORES_DIR", os.path.join(AIGCT_HOME, "processed_data"))

OUTPUT_DIR = os.environ.get("AIGCT_OUTPUT_DIR", os.path.join(AIGCT_HOME, "output"))

LABEL_INPUT = os.environ.get("AIGCT_LABEL_INPUT", os.path.join(OUTPUT_DIR, "input.txt"))
