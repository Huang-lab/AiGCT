# AIGCT benchmark curation

Code and inputs used to build the benchmark datasets distributed with
[AIGCT](https://github.com/Huang-lab/AiGCT), and the analysis notebooks that
produce Figures 2–4, Supplementary Figures S2–S5 and Supplementary Tables
S3–S4 of the AIGCT manuscript.

Figure 5 (hereditary cancer predisposition) is not reproduced here: it is
computed on individual-level UK Biobank exome and ICD-10 data, which is
controlled-access and cannot be redistributed. See the manuscript's Data
Availability statement.

This repository covers **how the benchmark database was made**. To *use* the
benchmark, install the `aigct` package instead — the curated database is
downloaded automatically and nothing here needs to be run.

```
curation/     pipeline that turns cohort variant lists into the benchmark tables
analysis/     scripts and notebooks that evaluate VEPs and produce the figures
data/         variant coordinate lists, transcript reference tables, gene lists
results/      supplementary tables and metric exports as published
config.yaml   paths to inputs and outputs
datasets.yaml the datasets that were built, and what each becomes in AIGCT
```

## What the pipeline does

Every benchmark task starts from a list of variant coordinates extracted from
the supplementary tables of published cohort studies (`data/annotation`). The
pipeline then:

1. **Annotates** each variant against dbNSFP v5.0a, pulling the rank scores of
   all 37 evaluated VEPs plus transcript, protein and gnomAD frequency fields
   (`curation/dbnsfp.py`).
2. **Resolves duplicate transcripts.** dbNSFP reports one record per
   transcript, so a variant in several transcripts appears several times.
   `curation/transcript_select.py` collapses each set to one representative
   transcript: the dbNSFP canonical transcript, else the CCDS transcript with
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
variants per gene to produce the class-balanced benchmark.

## External data you must download

Only the small reference tables are redistributed here. The large third-party
downloads are not, and must be fetched separately, then pointed at from
`config.yaml`:

| What | Where | Used by |
|---|---|---|
| dbNSFP v5.0a, academic release, unpacked (`dbNSFP5.0a_variant.chr*.gz`) | https://sites.google.com/site/jpopgen/dbNSFP | every stage |

Bundled under `data/reference/` for convenience:

| File | Source |
|---|---|
| `mart_export.txt` | Ensembl BioMart — transcript stable ID, version, length including UTRs and CDS, Ensembl canonical flag, CCDS ID |
| `CCDSID_length_table.current.csv` | CDS lengths derived from the NCBI CCDS release (`CCDS_nucleotide.current.fna.gz`) |

## Running it

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

## Analysis

`analysis/` queries the *published* AIGCT database rather than rebuilding it,
so it only needs the `aigct` package and its downloaded database.

- `generate_supp_tables.py` — Supplementary Tables S3 (AUC-ROC) and S4 (MWU),
  at both the 80% and 90% VEP-coverage thresholds. Set `AIGCT_CONFIG` to your
  `aigct.yaml` and `AIGCT_CLINVAR_CSV` to the ClinVar table from the curation
  step. S3 covers all 14 task–dataset combinations (28 sheets); S4 covers the
  10 non-ClinVar ones (20 sheets). MWU is not reported for ClinVar: n there is
  large enough that `-log10(p)` runs into the thousands — past the float64
  floor on the biggest strata — so it tracks sample size rather than effect
  size and ranks the VEPs no differently from AUC-ROC. ClinVar is assessed by
  AUC-ROC.
- `plot_vep.py` — the horizontal bar charts in Figures 2–4, with VEP labels
  coloured by training-data category (clinical-trained, population-tuned,
  population-free).
- `notebooks/` — per-task evaluation notebooks. Outputs are stripped before
  commit, so run a notebook top to bottom to reproduce its figures.

Nothing under `analysis/` hard-codes a filesystem location. The notebooks read
theirs from `notebooks/paths.py`, which takes every location from an
environment variable:

| Variable | Meaning | Default |
| --- | --- | --- |
| `AIGCT_HOME` | AIGCT installation holding `config/` and `db/` | this repository |
| `AIGCT_CONFIG` | `aigct.yaml` of that installation | `$AIGCT_HOME/config/aigct.yaml` |
| `AIGCT_DB_DIR` | benchmark database tables | `$AIGCT_HOME/db/data` |
| `AIGCT_CLINVAR_CSV` | ClinVar table from the curation step | `$AIGCT_HOME/output/processed/clinvar_withoutX.csv` |
| `AIGCT_SCORES_DIR` | MAVEN / EVE score tables | `$AIGCT_HOME/processed_data` |
| `AIGCT_OUTPUT_DIR` | where the notebooks write output | `$AIGCT_HOME/output` |

Setting `AIGCT_HOME` alone is usually enough:

```bash
export AIGCT_HOME=/path/to/your/aigct-install
```

`results/` holds the supplementary tables and metric exports exactly as
published.

## Notes on the data files

**Assemblies.** Source studies report on hg18, hg19 or hg38. The assembly of
each coordinate list is recorded in `datasets.yaml` and used to pick which
dbNSFP coordinate columns to join against.

**DDD gene lists.** The PrimateAI and AlphaMissense gene lists used to filter
the DDD task are not stored in this repository. They ship with the benchmark
database as `db/data/DDD/variant_filter_gene.csv`, under the filter codes
`DDD_RELATED_GENES_PRIMATEAI` (605 genes, from Sundaram et al.) and
`DDD_RELATED_GENES_ALPHA` (215 genes).

**Shared negative set.** ASD, CHD and DDD share one negative set, pooled from
the unaffected-sibling de novo variants of the ASD studies. Overlap removal
runs once per task, which is why the same five control files yield slightly
different counts for the three tasks in the released database.

**Superseded annotation files.** `data/annotation` retains earlier iterations
of several coordinate lists that are not referenced by `datasets.yaml`. They
are kept for provenance but were not used to build the released database:

- `ASD_case_annotation.txt`, `ASD_control_annotation.txt` — pooled ASD lists,
  replaced by the per-study `ASD_study{1..4}_*` files.
- `CHD_control_annotation.txt`, `DDD_control_annotation.txt`, `CHD_AIGCT.txt` —
  replaced by the shared control set.
- `MSK_passenger_annotation.txt` — the full MSK passenger set, narrowed to the
  subset in `MSK_passenger_hg19_annotation_6246.txt`.
- `TCGA_passenger_annotation.txt`, `TCGA_passenger_hg19_annotation_{1,2}.txt`,
  `TCGA_passenger_hg19_annotation_rd6000.txt` — replaced by the 5,000-variant
  random subset `TCGA_passenger_hg19_annotation_rd5000.txt`.
- `hotspot_annotation_fordb.txt`, `hotspot_annotation_fordb_withouttranscript.txt`
  — replaced by the driver-gene-filtered `..._gd.txt`.
- `DDD_alpha_case.txt`, `DDD_alpha_control.txt`, `ASD_primate_study2.txt` —
  intermediate gene-list-filtered subsets.

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
