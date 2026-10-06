"""Gene-level class balancing for the ClinVar benchmark.

Within each gene an equal number of pathogenic and benign variants is sampled,
using the smaller of the two groups as the target size, so that genes with a
lopsided ClinVar record cannot dominate the benchmark. Genes lacking either
label are excluded outright. Sampling within a gene is uniformly random under
a fixed seed.

The sampled set is exactly balanced (21,840 pathogenic and 21,840 benign
variants across 3,062 genes in the released database). The subset actually
evaluated is smaller and imbalanced, because the analyzer retains a variant
only if every VEP passing the coverage threshold scored it, and benign
variants are less completely covered by several VEPs.

Note on the released database. An earlier form of this module tried to draw
variants carrying an AlphaMissense score before those without one. That
preference never took effect: dbNSFP encodes a missing score as the string
"." and the check used `notna()`, which treats "." as present, so every
variant fell into the "scored" pool and the draw was a plain random sample.
The released set was built that way, and this module now does the plain
random sample explicitly; re-running it reproduces the released
`balanced_clinvar.csv` exactly.
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
    random_state=None,
):
    """Sample an equal number of positive and negative variants per gene."""
    cfg = load_config()["clinvar"]["balance"]
    gene_col = gene_col or cfg["gene_column"]
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

        balanced_group = pd.concat([
            pos.sample(n, random_state=random_state),
            neg.sample(n, random_state=random_state),
        ])
        balanced_data.append(balanced_group)
        gene_stats.append({"geneid": gene, "num_variants": len(balanced_group)})

    balanced_df = pd.concat(balanced_data).reset_index(drop=True)

    print("balanced variants:", len(balanced_df))
    print("genes:", len(gene_stats), " excluded (one label only):", len(excluded_genes))

    if gene_stats_path:
        pd.DataFrame(gene_stats).to_csv(gene_stats_path, index=False)
    if excluded_genes_path:
        pd.DataFrame({"excluded_geneid": excluded_genes}).to_csv(
            excluded_genes_path, index=False
        )
    balanced_df.to_csv(output, index=False)
    return balanced_df
