"""Generate every panel of Figures 2-4 and Supplementary Figures S2-S5.

Each panel is written as its own PNG (the multi-panel layouts in the manuscript
were assembled from these). Panel titles carry the evaluated counts in the same
form as the manuscript, e.g. ``Clinvar (balanced, +14704/-9544)``.

    AIGCT_CONFIG=/path/to/aigct.yaml python make_figures.py [--out DIR]
                                                            [--only NAME ...]
                                                            [--thresholds 80 90]

80% is the threshold used for Figures 2-4, 90% for Supplementary Figures S2-S4;
both are produced by default. Supplementary Figure S5 restricts the VEP panel to
the set benchmarked in the AlphaMissense study and is produced at 80% only.

Nothing beyond the `aigct` package and its downloaded benchmark database is
required: the ClinVar-excluded panels (Figure 4B/D/F) take their exclusion set
from the CLINVAR task's own label table.
"""
import argparse
import os
import sys
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_vep import VEP_TYPE, COLOR_MAP  # noqa: E402
from aigct.container import VEBenchmarkContainer  # noqa: E402
from aigct.model import VEQueryCriteria  # noqa: E402

warnings.filterwarnings("ignore")

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CONFIG = os.environ.get("AIGCT_CONFIG", os.path.join(REPO_ROOT, "config", "aigct.yaml"))
VARIANT_PK = ["GENOME_ASSEMBLY", "CHROMOSOME", "POSITION",
              "REFERENCE_NUCLEOTIDE", "ALTERNATE_NUCLEOTIDE"]

# MAVEN and MAVEN_(average) are present in the database but excluded from every
# analysis reported in the manuscript.
EXCLUDED_SOURCES = ["MAVEN", "MAVENAVG"]

# The VEP panels benchmarked in the AlphaMissense study, used for the
# cross-study validation in Supplementary Figure S5.
ALPHAMISSENSE_CLINVAR_PANEL = ["ALPHAM", "VAR_RL", "REVEL", "GMVP", "EIGEN", "CADD",
                               "POLYP2HVAR", "ESM1B", "SIFT", "POLYP2HDIV", "PRIMAI"]
ALPHAMISSENSE_CANCER_PANEL = ALPHAMISSENSE_CLINVAR_PANEL

# name -> (task, title stem, metric, query spec)
# "clinvar_excl" is resolved at run time to the CLINVAR label set.
PANELS = {
    # Figure 2 (80%) / Supplementary Figure S2 (90%) — AUC-ROC
    "Fig2A_ClinVar_balanced":     ("CLINVAR", "Clinvar (balanced",     "rocauc", dict(filter_names=["balanced_clinvar"])),
    "Fig2B_ClinVar_onestar":      ("CLINVAR", "Clinvar (Onestar+",     "rocauc", dict(filter_names=["Onestar", "Twostar", "Threestar", "Fourstar"])),
    "Fig2C_ClinVar_twostar":      ("CLINVAR", "Clinvar (Twostar+",     "rocauc", dict(filter_names=["Twostar", "Threestar", "Fourstar"])),
    "Fig2D_ClinVar_threestar":    ("CLINVAR", "Clinvar (Threestar+",   "rocauc", dict(filter_names=["Threestar", "Fourstar"])),
    # Figure 3 (80%) / Supplementary Figure S3 (90%) — AUC-ROC
    "Fig3A_Cancer_AlphaMissense": ("CANCER",  "Cancer (alphamissense", "rocauc", dict(filter_names=["Alpha_pos", "Alpha_neg"])),
    "Fig3B_Cancer_MSK":           ("CANCER",  "Cancer (MSK",           "rocauc", dict(filter_names=["MSK_hotspot", "MSK_passenger"])),
    "Fig3C_Cancer_TCGA":          ("CANCER",  "Cancer (TCGA",          "rocauc", dict(filter_names=["MSK_hotspot", "TCGA_passenger"])),
    "Fig3D_ADRD":                 ("ADRD",    "ADRD (",                "rocauc", None),
    # Figure 4 (80%) / Supplementary Figure S4 (90%) — MWU
    "Fig4A_ASD":                  ("ASD", "ASD (",                     "mwu", None),
    "Fig4B_ASD_exclClinVar":      ("ASD", "ASD (without Clinvar",      "mwu", "clinvar_excl"),
    "Fig4C_CHD":                  ("CHD", "CHD (",                     "mwu", None),
    "Fig4D_CHD_exclClinVar":      ("CHD", "CHD (without Clinvar",      "mwu", "clinvar_excl"),
    "Fig4E_DDD":                  ("DDD", "DDD (",                     "mwu", None),
    "Fig4F_DDD_exclClinVar":      ("DDD", "DDD (without Clinvar",      "mwu", "clinvar_excl"),
}

# Supplementary Figure S5 — fixed VEP panel, 80% threshold only.
S5_PANELS = {
    "FigS5A_ClinVar_balanced": ("CLINVAR", "balanced Clinvar (",
                                dict(filter_names=["balanced_clinvar"]), ALPHAMISSENSE_CLINVAR_PANEL),
    "FigS5B_Cancer_AlphaMissense": ("CANCER", "Cancer (AlphaMissense, ",
                                    dict(filter_names=["Alpha_pos", "Alpha_neg"]), ALPHAMISSENSE_CANCER_PANEL),
}


def draw(df, x_col, xlabel, title, path):
    """Horizontal bar chart with VEP labels coloured by training-data category."""
    df = df.sort_values(x_col, ascending=True).reset_index(drop=True)
    df["Category"] = df["SOURCE_NAME"].map(lambda v: VEP_TYPE.get(v.lower(), "Unknown"))
    fig, ax = plt.subplots(figsize=(8, 10))
    ax.barh(df["SOURCE_NAME"], df[x_col], color="#1f77b4", edgecolor="none", zorder=2)
    for label in ax.get_yticklabels():
        cat = df.loc[df["SOURCE_NAME"] == label.get_text(), "Category"]
        fc = COLOR_MAP.get(cat.iloc[0], "#FFFFFF") if len(cat) else "#FFFFFF"
        label.set_bbox({"facecolor": fc, "edgecolor": "none", "pad": 0})
        label.set_color("black")
        label.set_fontweight("bold")
        label.set_fontsize(18)
    if x_col == "ROC_AUC":
        ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
        ax.set_xlim(0, 1.1)
    ax.grid(axis="x", linestyle="--", color="gray", alpha=0.5, zorder=1)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xlabel(xlabel)
    ax.set_title(title, fontsize=20, loc="left")
    plt.tight_layout()
    plt.savefig(path, dpi=300)
    plt.close(fig)


def metric_frame(metrics, metric):
    if metric == "rocauc":
        return metrics.roc_metrics[["SOURCE_NAME", "ROC_AUC"]].dropna(), "ROC_AUC", "ROC AUC"
    return (metrics.mwu_metrics[["SOURCE_NAME", "NEG_LOG10_MWU_PVAL"]].dropna(),
            "NEG_LOG10_MWU_PVAL", "MWU -log10(pval)")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=os.path.join(REPO_ROOT, "results", "figures"))
    ap.add_argument("--only", nargs="*", help="panel names to produce (default: all)")
    ap.add_argument("--thresholds", nargs="*", type=int, default=[80, 90],
                    help="VEP coverage thresholds (80 = Figures 2-4, 90 = S2-S4)")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    container = VEBenchmarkContainer(CONFIG)
    clinvar_ids = container.query_mgr.get_variants_by_task("CLINVAR")[VARIANT_PK]
    clinvar_excl = VEQueryCriteria(variant_ids=clinvar_ids, include_variant_ids=False)

    rows = []
    for name, (task, stem, metric, q) in PANELS.items():
        if args.only and name not in args.only:
            continue
        qry = clinvar_excl if q == "clinvar_excl" else (VEQueryCriteria(**q) if q else None)
        for pct in args.thresholds:
            kw = dict(vep_min_overlap_percent=pct, variant_vep_retention_percent=100,
                      variant_effect_sources=EXCLUDED_SOURCES,
                      include_variant_effect_sources=False)
            if qry is not None:
                kw["variant_query_criteria"] = qry
            m = container.analyzer.compute_metrics(task, **kw)
            g = m.general_metrics.iloc[0]
            n_pos, n_neg = int(g.NUM_POSITIVE_LABELS), int(g.NUM_NEGATIVE_LABELS)
            sep = "" if stem.endswith("(") else ", "
            title = f"{stem}{sep}+{n_pos}/-{n_neg})"
            df, col, xlabel = metric_frame(m, metric)
            path = os.path.join(args.out, f"{name}_{pct}pct.png")
            draw(df, col, xlabel, title, path)
            rows.append(dict(panel=name, threshold=pct, n_pos=n_pos, n_neg=n_neg,
                             n_vep=len(df), file=os.path.basename(path)))
            print(f"{name} ({pct}%): {title}  {len(df)} VEPs")

    for name, (task, stem, q, panel) in S5_PANELS.items():
        if args.only and name not in args.only:
            continue
        m = container.analyzer.compute_metrics(
            task, variant_query_criteria=VEQueryCriteria(**q),
            vep_min_overlap_percent=80, variant_vep_retention_percent=100,
            variant_effect_sources=panel, include_variant_effect_sources=True)
        g = m.general_metrics.iloc[0]
        n_pos, n_neg = int(g.NUM_POSITIVE_LABELS), int(g.NUM_NEGATIVE_LABELS)
        title = f"{stem}+{n_pos:,} / -{n_neg:,})"
        df, col, xlabel = metric_frame(m, "rocauc")
        path = os.path.join(args.out, f"{name}.png")
        draw(df, col, xlabel, title, path)
        rows.append(dict(panel=name, threshold=80, n_pos=n_pos, n_neg=n_neg,
                         n_vep=len(df), file=os.path.basename(path)))
        print(f"{name}: {title}  {len(df)} VEPs")

    index = os.path.join(args.out, "panel_index.csv")
    pd.DataFrame(rows).to_csv(index, index=False)
    print(f"\nWrote {len(rows)} panels to {args.out}")


if __name__ == "__main__":
    main()
