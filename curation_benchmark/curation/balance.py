"""Gene-level class balancing for the ClinVar benchmark.

Within each gene an equal number of pathogenic and benign variants is sampled,
using the smaller of the two groups as the target size, so that genes with a
lopsided ClinVar record cannot dominate the benchmark. Genes lacking either
label are excluded outright.

Because the number of variants carrying a score from any given predictor
varies, variants that already have a score from a chosen priority VEP are
sampled first. Balancing happens before dbNSFP annotation is complete for
every sampled variant, so the resulting set is only approximately balanced —
see the Results section.
"""
import pandas as pd

from .config import load_config

LABEL_COLUMN = "BINARY_LABEL"


def balance_by_gene(
    df,
    output,
    gene_stats_path=None,
    excluded_genes_path=None,
    label_col=LABEL_COLUMN,
    gene_col=None,
    priority_col=None,
    random_state=None,
):
    """Sample an equal number of positive and negative variants per gene."""
    cfg = load_config()["clinvar"]["balance"]
    gene_col = gene_col or cfg["gene_column"]
    priority_col = priority_col or cfg["priority_column"]
    random_state = cfg["random_state"] if random_state is None else random_state

    # Only one- to four-star records are eligible.
    df = df[df["clinvar_review"].isin([1, 2, 3, 4])]

    balanced_data = []
    excluded_genes = []
    gene_stats = []

    for gene, group in df.groupby(gene_col):
        pos = group[group[label_col] == 1]
        neg = group[group[label_col] == 0]

        n = min(len(pos), len(neg))
        if n == 0:
            excluded_genes.append(gene)
            continue

        pos_with = pos[pos[priority_col].notna()]
        pos_without = pos[pos[priority_col].isna()]
        neg_with = neg[neg[priority_col].notna()]
        neg_without = neg[neg[priority_col].isna()]

        pos_sample = pd.concat([
            pos_with.sample(min(n, len(pos_with)), random_state=random_state),
            pos_without.sample(max(0, n - len(pos_with)), random_state=random_state),
        ]).head(n)

        neg_sample = pd.concat([
            neg_with.sample(min(n, len(neg_with)), random_state=random_state),
            neg_without.sample(max(0, n - len(neg_with)), random_state=random_state),
        ]).head(n)

        balanced_group = pd.concat([pos_sample, neg_sample])
        balanced_data.append(balanced_group)
        gene_stats.append({"geneid": gene, "num_variants": len(balanced_group)})

    balanced_df = pd.concat(balanced_data).reset_index(drop=True)
    balanced_df[priority_col] = balanced_df[priority_col].replace(pd.NA, ".")

    print("balanced variants:", len(balanced_df))
    print(
        f"without a {priority_col} score:",
        len(balanced_df[balanced_df[priority_col] == "."]),
    )

    if gene_stats_path:
        pd.DataFrame(gene_stats).to_csv(gene_stats_path, index=False)
    if excluded_genes_path:
        pd.DataFrame({"excluded_geneid": excluded_genes}).to_csv(
            excluded_genes_path, index=False
        )
    balanced_df.to_csv(output, index=False)
    return balanced_df
