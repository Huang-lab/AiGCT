"""
Creates db/aigct.db and populates three tables:
  - variant_effect_source   (from CSV)
  - variant_effect_task_auc (ROC AUC per task/VEP)
  - variant_effect_gene_auc (ROC AUC per task/VEP/gene)

Run from the repo root:
    python scripts/populate_db.py
"""

import os
import sys

# Allow imports from the repo root regardless of where the script is invoked.
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO_ROOT)

import pandas as pd
from sqlalchemy import create_engine, event
from sqlalchemy.orm import Session
from aigct.container import VEBenchmarkContainer
from aigct.db_model import (
    Base, VariantEffectSource, VariantEffectTaskAuc,
    VariantEffectGeneAuc)

VES_CSV = "/home/claudiof/gitrepo/aigct_data1/repo1/data/variant_effect_source.csv"
CONFIG_PATH = os.path.join(REPO_ROOT, "config", "aigct.yaml")
DB_PATH = os.path.join(REPO_ROOT, "db", "aigct.db")

TASK_CODES = ["ADRD", "CANCER", "CHD", "DDD", "ASD", "CLINVAR"]
# TASK_CODES = ["CLINVAR"]


def _enable_fk(dbapi_conn, _):
    dbapi_conn.execute("PRAGMA foreign_keys = ON")


def _load_variant_effect_sources(session: Session) -> None:
    df = pd.read_csv(VES_CSV)
    for _, row in df.iterrows():
        session.merge(
            VariantEffectSource(
                code=row["CODE"],
                name=row["NAME"],
                source_type=row["SOURCE_TYPE"],
                description=row["DESCRIPTION"],
            )
        )
    session.commit()
    print(f"  Loaded {len(df)} rows into variant_effect_source")


def _load_task_auc(session: Session, task_code: str, metrics) -> None:
    if metrics.roc_metrics is None or metrics.general_metrics is None:
        print(f"  [{task_code}] No ROC or general metrics — skipping task AUC")
        return

    df = metrics.roc_metrics.merge(
        metrics.general_metrics[["SCORE_SOURCE", "NUM_POSITIVE_LABELS", "NUM_NEGATIVE_LABELS"]],
        on="SCORE_SOURCE",
        how="left",
    )
    df = df.merge(
        metrics.mwu_metrics[["SCORE_SOURCE", "NEG_LOG10_MWU_PVAL"]],
        on="SCORE_SOURCE",
        how="left",
    )
    for _, row in df.iterrows():
        session.merge(
            VariantEffectTaskAuc(
                task_code=task_code,
                score_source=row["SCORE_SOURCE"],
                auc=row["ROC_AUC"] if pd.notna(row["ROC_AUC"]) else None,
                neg_log10_mwu_pval=row["NEG_LOG10_MWU_PVAL"] if pd.notna(row["NEG_LOG10_MWU_PVAL"]) else None,
                num_positive=int(row["NUM_POSITIVE_LABELS"]),
                num_negative=int(row["NUM_NEGATIVE_LABELS"]),
            )
        )
    session.commit()
    print(f"  [{task_code}] Loaded {len(df)} rows into variant_effect_task_auc")


def _load_gene_auc(session: Session, task_code: str, metrics) -> None:
    if metrics.gene_roc_metrics is None or metrics.gene_general_metrics is None:
        print(f"  [{task_code}] No gene ROC or general metrics — skipping gene AUC")
        return

    df = metrics.gene_roc_metrics.merge(
        metrics.gene_general_metrics[
            ["SCORE_SOURCE", "GENE_SYMBOL", "NUM_POSITIVE_LABELS", "NUM_NEGATIVE_LABELS"]
        ],
        on=["SCORE_SOURCE", "GENE_SYMBOL"],
        how="left",
    )
    df = df.merge(
        metrics.gene_mwu_metrics[["SCORE_SOURCE", "GENE_SYMBOL", "NEG_LOG10_MWU_PVAL"]],
        on=["SCORE_SOURCE", "GENE_SYMBOL"],
        how="left",
    )
    for _, row in df.iterrows():
        session.merge(
            VariantEffectGeneAuc(
                task_code=task_code,
                score_source=row["SCORE_SOURCE"],
                gene_symbol=row["GENE_SYMBOL"],
                auc=row["ROC_AUC"] if pd.notna(row["ROC_AUC"]) else None,
                neg_log10_mwu_pval=row["NEG_LOG10_MWU_PVAL"] if pd.notna(row["NEG_LOG10_MWU_PVAL"]) else None,
                num_positive=int(row["NUM_POSITIVE_LABELS"]),
                num_negative=int(row["NUM_NEGATIVE_LABELS"]),
            )
        )
    session.commit()
    print(f"  [{task_code}] Loaded {len(df)} rows into variant_effect_gene_auc")


def main():
    engine = create_engine(f"sqlite:///{DB_PATH}")
    event.listen(engine, "connect", _enable_fk)
    Base.metadata.create_all(engine)
    print(f"Database created at {DB_PATH}")

    container = VEBenchmarkContainer(CONFIG_PATH)
    analyzer = container.analyzer

    with Session(engine) as session:
        print("\nDeleting existing rows...")
        session.query(VariantEffectGeneAuc).delete()
        session.query(VariantEffectTaskAuc).delete()
        session.query(VariantEffectSource).delete()
        session.commit()

        print("\nLoading variant_effect_source...")
        _load_variant_effect_sources(session)

        for task_code in TASK_CODES:
            print(f"\nComputing metrics for {task_code}...")
            metrics = analyzer.compute_metrics(
                task_code,
                compute_gene_metrics=True,
                vep_min_overlap_percent=80,
                variant_vep_retention_percent=100,
                # metrics="mwu"
            )
            _load_task_auc(session, task_code, metrics)
            _load_gene_auc(session, task_code, metrics)

    print("\nDone.")


if __name__ == "__main__":
    main()
