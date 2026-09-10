import matplotlib.pyplot as plt

VEP_TYPE = {
    "alphamissense":        "population-tuned",
    "bayesdel_addaf":       "clinical-trained",
    "bayesdel_noaf":        "clinical-trained",
    "cadd_raw":             "clinical-trained",
    "clinpred":             "clinical-trained",
    "dann":                 "clinical-trained",
    "deogen2":              "clinical-trained",
    "eigen-raw_coding":     "clinical-trained",
    "eigen-pc-raw_coding":  "clinical-trained",
    "esm1b":                "population-free",
    "fathmm-xf_coding":     "clinical-trained",
    "gmvp":                 "clinical-trained",
    "list-s2":              "population-tuned",
    "m-cap":                "clinical-trained",
    "maven":                "population-free",
    "maven_(average)":      "population-free",
    "metalr":               "clinical-trained",
    "metarnn":              "clinical-trained",
    "metasvm":              "clinical-trained",
    "mpc":                  "clinical-trained",
    "mutationassessor":     "population-free",
    "mutationtaster":       "clinical-trained",
    "mutformer":            "clinical-trained",
    "mutpred":              "clinical-trained",
    "mutscore":             "clinical-trained",
    "mvp":                  "clinical-trained",
    "phactboost":           "clinical-trained",
    "polyphen2_hdiv":       "clinical-trained",
    "polyphen2_hvar":       "clinical-trained",
    "primateai":            "population-tuned",
    "provean":              "population-free",
    "revel":                "clinical-trained",
    "sift":                 "population-free",
    "sift4g":               "population-free",
    "varity_er":            "clinical-trained",
    "varity_er_loo":        "clinical-trained",
    "varity_r":             "clinical-trained",
    "varity_r_loo":         "clinical-trained",
    "vest4":                "clinical-trained",
}

COLOR_MAP = {
    "population-tuned": "#D62728",  # red
    "clinical-trained": "#FFD700",  # yellow
    "population-free":  "#2CA02C",  # green
    "Unknown":          "#808080",
}


def plot_vep_results(metrics, metric="rocauc", output_path="vep_auc_colored.png"):
    """
    Parameters
    ----------
    metrics    : VEAnalysisResult from container.analyzer.compute_metrics()
    metric     : "rocauc" or "mwu"
    output_path: where to save the figure
    """
    if metric == "rocauc":
        df = metrics.roc_metrics[["SOURCE_NAME", "ROC_AUC"]].dropna(subset=["ROC_AUC"]).copy()
        x_col  = "ROC_AUC"
        xlabel = "ROC AUC"
    elif metric == "mwu":
        df = metrics.mwu_metrics[["SOURCE_NAME", "NEG_LOG10_MWU_PVAL"]].dropna(subset=["NEG_LOG10_MWU_PVAL"]).copy()
        x_col  = "NEG_LOG10_MWU_PVAL"
        xlabel = "MWU -log10(pval)"
    else:
        raise ValueError(f"metric must be 'rocauc' or 'mwu', got '{metric}'")

    df = df.sort_values(x_col, ascending=True).reset_index(drop=True)
    df["Category"] = df["SOURCE_NAME"].apply(
        lambda v: next((VEP_TYPE[k] for k in VEP_TYPE if k == v.lower()), "Unknown")
    )

    fig, ax = plt.subplots(figsize=(8, 10))

    # bars: uniform blue
    ax.barh(df["SOURCE_NAME"], df[x_col], color="#1f77b4", edgecolor="none", zorder=2)

    # label backgrounds: colored by VEP category
    for label in ax.get_yticklabels():
        rows = df.loc[df["SOURCE_NAME"] == label.get_text(), "Category"].values
        fc = COLOR_MAP.get(rows[0], "#FFFFFF") if len(rows) else "#FFFFFF"
        label.set_bbox({"facecolor": fc, "edgecolor": "none", "pad": 0})
        label.set_color("black")
        label.set_fontweight("bold")
        label.set_fontsize(18)

    if metric == "rocauc":
        ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
        ax.set_xlim(0, 1.1)

    ax.grid(axis="x", linestyle="--", color="gray", alpha=0.5, zorder=1)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xlabel(xlabel)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.show()
    print(f"Saved → {output_path}")
