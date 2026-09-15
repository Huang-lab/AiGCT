"""
Generate supplementary tables from AIGCT benchmark results.
Runs both parameter sets:
  - 80%  overlap  (vep_min_overlap_percent=80)
  - 90%  overlap  (vep_min_overlap_percent=90)
Both use variant_vep_retention_percent=100.

Output:
  supp_table_ROC_AUC.xlsx  — Supplementary Table S3, one sheet per task per
                             threshold (14 task-dataset combinations x 2)
  supp_table_MWU.xlsx      — Supplementary Table S4, same layout but without
                             the four ClinVar strata (10 combinations x 2);
                             see the MWU block below for why
  supp_table_PR_AUC.xlsx   — auPRC companion table, same layout as S3
  supp_table_dataset_counts.xlsx — per task: database totals, the counts that
                             survive the coverage/retention filters, and how
                             many VEPs passed the threshold

Covers Figures 2, 3 and 4. The hereditary cancer predisposition analysis
(Figure 5) uses AUBPRC on controlled-access UK Biobank data and is not
generated here.
"""
import pandas as pd
import os
from aigct.container import VEBenchmarkContainer
from aigct.model import VEQueryCriteria

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Point AIGCT_CONFIG at the aigct.yaml of an AIGCT installation that has already
# downloaded the benchmark database. Nothing else is required.
CONFIG = os.environ.get("AIGCT_CONFIG", os.path.join(REPO_ROOT, "config", "aigct.yaml"))
OUTPUT_DIR = os.environ.get(
    "AIGCT_SUPP_TABLE_DIR", os.path.join(REPO_ROOT, "results", "supp_tables")
)
os.makedirs(OUTPUT_DIR, exist_ok=True)

container = VEBenchmarkContainer(CONFIG)

# Figures 4B/D/F exclude every variant carrying a ClinVar classification. That
# set is exactly the label table of the CLINVAR task, so it is read from the
# database rather than from a curation output file.
VARIANT_PK = ["GENOME_ASSEMBLY", "CHROMOSOME", "POSITION",
              "REFERENCE_NUCLEOTIDE", "ALTERNATE_NUCLEOTIDE"]
clinvar_df = container.query_mgr.get_variants_by_task("CLINVAR")[VARIANT_PK]
clinvar_excl = VEQueryCriteria(variant_ids=clinvar_df, include_variant_ids=False)

# ── Task definitions ──────────────────────────────────────────────────────────
# Each entry: (sheet_name, task_code, VEQueryCriteria_or_None)
TASKS = [
    ("ADRD",                  "ADRD",    None),
    ("ASD",                   "ASD",     None),
    ("ASD_excl_ClinVar",      "ASD",     clinvar_excl),
    ("CHD",                   "CHD",     None),
    ("CHD_excl_ClinVar",      "CHD",     clinvar_excl),
    ("DDD",                   "DDD",     None),
    ("DDD_excl_ClinVar",      "DDD",     clinvar_excl),
    ("Cancer_MSK",            "CANCER",  VEQueryCriteria(filter_names=["MSK_hotspot", "MSK_passenger"])),
    ("Cancer_TCGA",           "CANCER",  VEQueryCriteria(filter_names=["MSK_hotspot", "TCGA_passenger"])),
    ("Cancer_AlphaMissense",  "CANCER",  VEQueryCriteria(filter_names=["Alpha_pos", "Alpha_neg"])),
    ("ClinVar_Balanced",      "CLINVAR", VEQueryCriteria(filter_names=["balanced_clinvar"])),
    # ClinVar review-status strata (Figure 2B-D). Each is cumulative: "Onestar+"
    # means one star or better, so it is the union of the one-, two-, three- and
    # four-star filters defined in data/CLINVAR/variant_filter.csv.
    ("ClinVar_Onestar+",      "CLINVAR", VEQueryCriteria(
        filter_names=["Onestar", "Twostar", "Threestar", "Fourstar"])),
    ("ClinVar_Twostar+",      "CLINVAR", VEQueryCriteria(
        filter_names=["Twostar", "Threestar", "Fourstar"])),
    ("ClinVar_Threestar+",    "CLINVAR", VEQueryCriteria(
        filter_names=["Threestar", "Fourstar"])),
]

THRESHOLDS = [80, 90]

# ── Run and collect results ───────────────────────────────────────────────────
roc_sheets = {}
mwu_sheets = {}
pr_sheets = {}
count_rows = []

for overlap in THRESHOLDS:
    for sheet_name, task_code, qry in TASKS:
        label = f"{sheet_name} ({overlap}%)"
        print(f"Running: {label} ...")

        kwargs = dict(
            vep_min_overlap_percent=overlap,
            variant_vep_retention_percent=100,
        )
        if qry is not None:
            kwargs["variant_query_criteria"] = qry

        kwargs["variant_effect_sources"] = ["MAVEN", "MAVENAVG"]
        kwargs["include_variant_effect_sources"] = False

        metrics = container.analyzer.compute_metrics(task_code, **kwargs)

        n      = int(metrics.general_metrics["NUM_VARIANTS"].iloc[0])
        n_pos  = int(metrics.general_metrics["NUM_POSITIVE_LABELS"].iloc[0])
        n_neg  = int(metrics.general_metrics["NUM_NEGATIVE_LABELS"].iloc[0])

        # ROC AUC
        roc_df = (
            metrics.roc_metrics[["SOURCE_NAME", "ROC_AUC"]]
            .rename(columns={"SOURCE_NAME": "VEP"})
            .sort_values("ROC_AUC", ascending=False)
            .reset_index(drop=True)
        )
        roc_df.insert(0, "Variant_Total",    n)
        roc_df.insert(1, "Positive_Labels",  n_pos)
        roc_df.insert(2, "Negative_Labels",  n_neg)
        roc_sheets[label] = roc_df

        # auPRC
        pr_df = (
            metrics.pr_metrics[["SOURCE_NAME", "PR_AUC"]]
            .rename(columns={"SOURCE_NAME": "VEP", "PR_AUC": "auPRC"})
            .sort_values("auPRC", ascending=False)
            .reset_index(drop=True)
        )
        pr_df.insert(0, "Variant_Total",    n)
        pr_df.insert(1, "Positive_Labels",  n_pos)
        pr_df.insert(2, "Negative_Labels",  n_neg)
        pr_sheets[label] = pr_df

        # MWU — not reported for ClinVar. On these strata n is large enough that
        # -log10(p) runs into the thousands and, for the biggest ones, past the
        # float64 floor; the statistic tracks sample size rather than effect
        # size and orders the VEPs no differently from AUC-ROC. ClinVar is
        # assessed by AUC-ROC (Figures 2 and S2); MWU is reported for the tasks
        # it was introduced for, where label noise attenuates rank-based
        # discrimination.
        if task_code != "CLINVAR":
            mwu_df = (
                metrics.mwu_metrics[["SOURCE_NAME", "NEG_LOG10_MWU_PVAL"]]
                .rename(columns={"SOURCE_NAME": "VEP",
                                 "NEG_LOG10_MWU_PVAL": "MWU_neg_log10_pval"})
                .sort_values("MWU_neg_log10_pval", ascending=False)
                .reset_index(drop=True)
            )
            mwu_df.insert(0, "Variant_Total",    n)
            mwu_df.insert(1, "Positive_Labels",  n_pos)
            mwu_df.insert(2, "Negative_Labels",  n_neg)
            mwu_sheets[label] = mwu_df

        # Database totals for the task, before any coverage filtering.
        task_variants = container.query_mgr.get_variants_by_task(task_code)
        count_rows.append({
            "Task": task_code,
            "Dataset": sheet_name,
            "Coverage_threshold_pct": overlap,
            "db_variants": len(task_variants),
            "db_positive": int((task_variants["BINARY_LABEL"] == 1).sum()),
            "db_negative": int((task_variants["BINARY_LABEL"] == 0).sum()),
            # Ensembl gene IDs, matching the gene counts in Figure 1B.
            "db_genes": task_variants["ENSEMBL_GENE_ID"].nunique(),
            "evaluated_variants": n,
            "evaluated_positive": n_pos,
            "evaluated_negative": n_neg,
            "VEPs_passing_threshold": len(roc_df),
        })

        print(f"  → {n} variants ({n_pos}+/{n_neg}-), {len(roc_df)} VEPs")

# ── Save Excel ────────────────────────────────────────────────────────────────
roc_path = os.path.join(OUTPUT_DIR, "supp_table_ROC_AUC.xlsx")
mwu_path = os.path.join(OUTPUT_DIR, "supp_table_MWU.xlsx")
pr_path  = os.path.join(OUTPUT_DIR, "supp_table_PR_AUC.xlsx")
cnt_path = os.path.join(OUTPUT_DIR, "supp_table_dataset_counts.xlsx")

with pd.ExcelWriter(roc_path, engine="openpyxl") as writer:
    for sheet_label, df in roc_sheets.items():
        df.to_excel(writer, sheet_name=sheet_label, index=False)

with pd.ExcelWriter(mwu_path, engine="openpyxl") as writer:
    for sheet_label, df in mwu_sheets.items():
        df.to_excel(writer, sheet_name=sheet_label, index=False)

with pd.ExcelWriter(pr_path, engine="openpyxl") as writer:
    for sheet_label, df in pr_sheets.items():
        df.to_excel(writer, sheet_name=sheet_label, index=False)

print(f"\nSaved: {roc_path}")
print(f"Saved: {mwu_path}")
with pd.ExcelWriter(cnt_path, engine="openpyxl") as writer:
    pd.DataFrame(count_rows).to_excel(writer, sheet_name="dataset_counts", index=False)

print(f"Saved: {pr_path}")
print(f"Saved: {cnt_path}")
print(f"ROC sheets ({len(roc_sheets)}): {list(roc_sheets.keys())}")
print(f"MWU sheets ({len(mwu_sheets)}): {list(mwu_sheets.keys())}")
