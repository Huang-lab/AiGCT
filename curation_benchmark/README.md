# AIGCT benchmark curation

Code and inputs used to build the benchmark datasets distributed with
[AIGCT](https://github.com/Huang-lab/AiGCT), and the analysis scripts that
produce Figures 2–4, Supplementary Figures S2–S5 and Supplementary Tables
S3–S4 of the AIGCT manuscript.

This directory contains only what is needed to go from dbNSFP to the published
results. Figure 5 (hereditary cancer predisposition) is not reproduced here: it
is computed on individual-level UK Biobank exome and ICD-10 data, which is
controlled-access and cannot be redistributed. See the manuscript's Data
Availability statement.

This covers **how the benchmark database was made**. To *use* the benchmark,
install the `aigct` package instead — the curated database is downloaded
automatically and nothing here needs to be run.

```
curation/     pipeline that turns cohort variant lists into the benchmark tables
analysis/     scripts that evaluate VEPs and produce the figures and tables
data/         variant coordinate lists, transcript reference tables, gene lists
results/      supplementary tables as published
config.yaml   paths to inputs and outputs
datasets.yaml the datasets that were built, and what each becomes in AIGCT
```

How the curated CSVs become the distributed database is documented under
"From curated CSVs to the benchmark database" below.

## What the pipeline does

Every benchmark task starts from a list of variant coordinates extracted from
the supplementary tables of published cohort studies (`data/annotation`). The
pipeline then:

1. **Annotates** each variant against dbNSFP v5.0a, pulling the rank scores of
   all 37 evaluated VEPs (Supplementary Table S2) plus transcript, protein and
   gnomAD frequency fields (`curation/dbnsfp.py`).
2. **Resolves duplicate transcripts.** dbNSFP reports one record per
   transcript, so a variant in several transcripts appears several times.
   `curation/transcript_select.py` collapses each set to one representative
   transcript: the dbNSFP canonical transcript where there is exactly one;
   otherwise the Ensembl canonical transcript, else the CCDS transcript with
   the longest CDS (ties broken by overall transcript length), else the longest
   transcript overall. This is the protocol in Supplementary Figure S6.
3. **Removes overlaps** between the positive and negative set of each task, so
   that a de novo variant seen in both a case cohort and the shared control set
   is dropped from both (`curation/overlap.py`).

The ClinVar task is different: it has no coordinate list. `curation/clinvar.py`
scans the whole dbNSFP release and keeps every record carrying a ClinVar
classification, using dbNSFP's own `clinvar_clnsig` and `clinvar_review`
fields — so **dbNSFP fixes which ClinVar release the benchmark reflects**.
`curation/balance.py` then samples an equal number of pathogenic and benign
variants per gene to produce the class-balanced benchmark (21,840 + 21,840
variants across 3,062 genes).

### VEP versions are fixed by the dbNSFP release

Because every VEP score comes from one dbNSFP release, the model version of
each predictor is whatever that release ships. dbNSFP v5.0a supplies **CADD
v1.7** and **MutationTaster2021**, among others. This matters when comparing
against previously published benchmarks: the AlphaMissense study (Cheng et al.,
2023) predates dbNSFP v4.7 and therefore used **CADD v1.6**, a model without
the protein language model and regulatory features added in v1.7. CADD's
position in the cross-study comparison (Supplementary Figure S5) reflects a
different predictor, not a difference in the evaluation.

### What changed between database release 1.0.0 and 1.0.1

Release 1.0.0 carried no CADD_raw, Eigen-raw_coding or Eigen-PC-raw_coding
scores for the CLINVAR task. The dbNSFP scan did extract them, but wrote the
raw score columns under their bare dbNSFP names (`CADD_raw`, …) while the
loader expected the renamed form (`CADD_raw_score`, …) and skipped what it
could not find. `curation/clinvar.py` applies the rename (`VEP_RENAME` in
`curation/columns.py`), so a fresh run is unaffected.

1.0.0 also carried MAVEN and MAVEN_(average) scores, which are not produced by
this pipeline and are not among the predictors evaluated in the manuscript.

Release 1.0.1 therefore adds the three missing predictors to the ClinVar task
and drops the two MAVEN sources, leaving exactly the **37 VEPs of
Supplementary Table S2** in every task. Dropping MAVEN changes no reported
number: it was excluded from every analysis in 1.0.0 as well, and the
supplementary tables regenerate bit for bit from either release.

## From curated CSVs to the benchmark database

The curation stages above end at one annotated, transcript-resolved CSV per
dataset under `output/processed/`. Those CSVs are then loaded into the
flat-file relational database that the `aigct` package reads. This repository
documents that step rather than shipping a driver for it, because the released
database is downloaded from Zenodo and nothing here needs to rebuild it.

The database is a directory of CSV files: four global tables plus one
subdirectory per task.

| Table | Scope | Contents |
|---|---|---|
| `variant.csv` | global | master list of variants, keyed by (assembly, chromosome, position, ref, alt), with hg19/hg18 coordinates, amino acid change, gene/transcript/protein IDs and gnomAD allele frequency |
| `variant_task.csv` | global | the six task codes |
| `variant_effect_source.csv` | global | the 37 VEP codes and display names |
| `variant_data_source.csv` | global | allele-frequency sources |
| `<TASK>/variant_effect_label.csv` | per task | one row per variant: `LABEL_SOURCE`, `BINARY_LABEL` |
| `<TASK>/variant_effect_score.csv` | per task | one row per variant **per VEP**: `SCORE_SOURCE`, `RAW_SCORE`, `RANK_SCORE` |
| `<TASK>/variant_filter.csv` | per task | the named filters available for that task |
| `<TASK>/variant_filter_variant.csv` | per task | which variants belong to each variant-based filter |
| `<TASK>/variant_filter_gene.csv` | per task | which genes belong to each gene-based filter |

All variants are stored on **hg38**; coordinate lists reported on hg19 or hg18
are carried in the `PRIOR_*` columns rather than as separate rows.

**Loading a curated CSV.** `aigct.etl.repo_loader.RepositoryLoader.load_variant_file`
takes one processed CSV together with the task, a `LABEL_SOURCE` and a binary
label, and expands it into three tables at once: one row in `variant.csv`, one
row in `variant_effect_label.csv`, and one row in `variant_effect_score.csv`
for every VEP score column present. Its `VEP_COLUMN_LIST` is what maps a
dbNSFP column pair such as `AlphaMissense_score` / `AlphaMissense_rankscore`
onto the database VEP code `ALPHAM`. Each call carries a fixed label, so
positive and negative sets are loaded separately — which is why
`datasets.yaml` lists case and control files as separate datasets.

`LABEL_SOURCE` records which curated set a variant came from. Five of the six
tasks use a single source named after the task; the cancer task distinguishes
its five contributing sets:

```
ADRD / ASD / CHD / DDD / CLINVAR : one source per task
CANCER : HOTSPOT, ALPHAMISSENSE_POS (positives)
         MSK_PASSENGER, TCGA_PASSENGER, ALPHAMISSENSE_NEG (negatives)
```

**Named filters.** The filters that Figures 2–4 select on are loaded
separately from the labels, so a variant can carry one label and belong to
several filters. Every `filter_code` in `datasets.yaml` becomes a row in the
task's `variant_filter.csv` and a set of rows in
`variant_filter_variant.csv`. Gene-based filters (the cancer driver-gene
categories and the DDD gene lists) live in `variant_filter_gene.csv` instead.
As released:

| Task | Filters |
|---|---|
| CANCER | `MSK_HOTSPOT` (826), `ALPHA_POS` (862), `MSK_PASSENGER` (6,234), `TCGA_PASSENGER` (5,000), `ALPHA_NEG` (1,733); gene filters `TSG` (99 genes), `ONCOG` (84 genes) |
| CLINVAR | `ONESTAR` (122,853), `TWOSTAR` (46,661), `THREESTAR` (3,161), `FOURSTAR` (8), `BALANCED_CLINVAR` (43,680) |
| ASD | `ASD_CASE1`–`4`, `ASD_PRI_CASE1`–`3`, `ASD_CONTROL1`–`4`, `ASD_PRI_CONTROL2` |
| CHD | `CHD_CASE`, `CHD_CONTROL1`–`4`, `CHD_PRI_CONTROL2` |
| DDD | `DDD_CASE`, `DDD_PRI_CASE1`–`2`, `DDD_CONTROL1`–`4`, `DDD_PRI_CONTROL2`; gene filters `DDD_RELATED_GENES_PRIMATEAI` (605 genes), `DDD_RELATED_GENES_ALPHA` (215 genes) |
| ADRD | none — the task is a single case/control set |

Two consequences of this layout are worth spelling out, because the
manuscript's numbers depend on them:

- **The ClinVar star filters hold one star rating each, not a cumulative
  range.** `ONESTAR` is the one-star variants only, and the four filters
  partition the task exactly (122,853 + 46,661 + 3,161 + 8 = 172,683). The
  "one star or higher" stratum in Figure 2B is therefore the *union* of all
  four filters, which is how `generate_supp_tables.py` and `make_figures.py`
  request it.
- **A variant can sit in more than one filter but carries only one label.**
  In the cancer task the 862 `ALPHA_POS` variants overlap the 826
  `MSK_HOTSPOT` variants by 675, so only 187 of them are labelled
  `ALPHAMISSENSE_POS`; the union is the 1,013 positives quoted in the
  Methods (826 + 862 − 675).

**Note on the loader.** `aigct.etl.repo_loader` in the published package
provides `init_variant_task`, `init_variant_effect_source` and
`load_variant_file`. The filter tables and the ClinVar task of the released
database were built with additional loader routines that are not part of the
published package, so the database cannot currently be rebuilt end to end from
PyPI alone. This does not affect use of the benchmark: the database itself is
distributed through Zenodo and installed by `install_db`.

## External data you must download

Only the small reference tables are redistributed here. The large third-party
download is not, and must be fetched separately, then pointed at from
`config.yaml`:

| What | Where | Used by |
|---|---|---|
| dbNSFP v5.0a, academic release, unpacked (`dbNSFP5.0a_variant.chr*.gz`) | https://sites.google.com/site/jpopgen/dbNSFP | every stage |

Bundled under `data/reference/` for convenience:

| File | Source |
|---|---|
| `mart_export.txt` | Ensembl BioMart — transcript stable ID, version, length including UTRs and CDS, Ensembl canonical flag, CCDS ID |
| `CCDSID_length_table.current.csv` | CDS lengths derived from the NCBI CCDS release (`CCDS_nucleotide.current.fna.gz`) |

## Running the curation

```bash
pip install -r requirements.txt

# edit config.yaml so external.dbnsfp_dir points at your dbNSFP release
python -m curation.run_curation annotate          # all cohort datasets
python -m curation.run_curation annotate ASD_case1  # or just one
python -m curation.run_curation clinvar           # scan dbNSFP for ClinVar
python -m curation.run_curation balance           # gene-balance ClinVar
python -m curation.run_curation overlap           # positive/negative split
python -m curation.run_curation exclude-clinvar   # Figure 4B/D/F inputs
python -m curation.run_curation all               # everything, in order
```

Every stage writes under `output/` (git-ignored) and is safe to re-run. Both
`annotate` and `clinvar` stream the full dbNSFP release from disk and take
hours; `clinvar` reads all 25 chromosome files end to end.

`datasets.yaml` is the manifest: it lists each dataset, the coordinate list it
comes from, the assembly that list is reported on, the CSV it produces, and the
`filter_code` it becomes in the AIGCT database. Every benchmark number in the
manuscript can be traced back through it.

### Tests

```bash
pytest curation_benchmark/tests
```

Unit tests for the two stages where a defect is silent rather than loud: the
transcript hierarchy, which must never return an empty group (a variant deleted
there leaves no message and no duplicates entry), and the gene balancing, whose
printed coverage count is the only readout of how well priority-first sampling
worked. They need neither dbNSFP nor the reference tables, so they run in
seconds on a bare checkout.

The repository root's `pytest.ini` sets `testpaths=tests`, which is the `aigct`
package's own suite and needs a downloaded benchmark database; these are kept
separate and named explicitly.

## Running the analysis

`analysis/` queries the *published* AIGCT database rather than rebuilding it,
so it needs only the `aigct` package and its downloaded database. Point
`AIGCT_CONFIG` at the `aigct.yaml` of that installation; nothing else is
required.

```bash
export AIGCT_CONFIG=/path/to/your/aigct-install/config/aigct.yaml

python analysis/generate_supp_tables.py   # Supplementary Tables S3, S4 and counts
python analysis/make_figures.py           # all panels of Figures 2–4, S2–S5
```

- `generate_supp_tables.py` — Supplementary Tables S3 (AUC-ROC) and S4 (MWU) at
  both the 80% and 90% VEP-coverage thresholds, plus the auPRC companion table
  and a dataset-counts table. S3 covers all 14 task–dataset combinations (28
  sheets); S4 covers the 10 non-ClinVar ones (20 sheets). MWU is not reported
  for ClinVar: n there is large enough that `-log10(p)` runs into the thousands
  — past the float64 floor on the biggest strata — so it tracks sample size
  rather than effect size and ranks the VEPs no differently from AUC-ROC.
  ClinVar is assessed by AUC-ROC.
- `make_figures.py` — every panel of Figures 2–4 (80% threshold) and
  Supplementary Figures S2–S4 (90%), plus the two cross-study panels of
  Supplementary Figure S5, written as individual PNGs with the evaluated counts
  in each panel title. `--only` restricts to named panels, `--thresholds` to one
  threshold.
- `plot_vep.py` — the bar-chart styling shared by those panels, with VEP labels
  coloured by training-data category (clinical-trained, population-tuned,
  population-free).

`results/supp_tables/` holds the supplementary tables exactly as published.

## Notes on the data files

**Assemblies.** Source studies report on hg18, hg19 or hg38. The assembly of
each coordinate list is recorded in `datasets.yaml` and used to pick which
dbNSFP coordinate columns to join against.

**DDD gene lists.** The PrimateAI and AlphaMissense gene lists associated with
the DDD task ship with the benchmark database as
`db/data/DDD/variant_filter_gene.csv`, under the filter codes
`DDD_RELATED_GENES_PRIMATEAI` (605 genes, from Sundaram et al.) and
`DDD_RELATED_GENES_ALPHA` (215 genes). They are available as named gene filters
but were **not applied** in any analysis reported in the manuscript, which uses
all DDD de novo missense variants.

**Shared negative set.** ASD, CHD and DDD share one negative set, pooled from
the unaffected-sibling de novo variants of the ASD studies. Overlap removal
runs once per task, which is why the same five control files yield slightly
different counts for the three tasks in the released database.

**TCGA passengers.** The full TCGA passenger set was randomly reduced: a draw
of 6,000 variants, of which 5,909 had a dbNSFP record, from which a random
5,000 were retained. `data/annotation/TCGA_passenger_hg19_annotation_rd5000.txt`
lists exactly those 5,000 (hg19), so the pipeline reproduces the published
`TCGA_PASSENGER` filter without re-sampling.

**One reconstructed file.** `MSK_passenger_hg19_annotation_6246.txt` was
regenerated from the `hg19_chr`, `hg19_pos(1-based)`, `ref` and `alt` columns
of the curated `MSK_passenger_6246.csv`, because the original coordinate list
had been destroyed: an earlier version of the annotation step took the input
path and the "variants not found in dbNSFP" output path as separate arguments,
and this dataset was run with the same path passed for both, truncating the
input to a bare header. The reconstruction holds 6,239 variants — those that
were found in dbNSFP — rather than the 6,246 the filename refers to. The seven
missing variants had no dbNSFP record and so never reached the released
database; the curated output, and every benchmark result derived from it, is
unaffected.
