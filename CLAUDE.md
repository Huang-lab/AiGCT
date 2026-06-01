# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**AIGCT (AI Genomics CollecTive)** is a Python platform for benchmarking variant effect predictors (VEPs) across genomics-based precision medicine tasks (CANCER, ADRD, CHD, DDD, ASD, CLINVAR). It compares VEP predictions against ground-truth variant labels and produces metrics (ROC-AUC, PR-AUC, Mann-Whitney U), visualizations, and CSV exports.

## Commands

```bash
# Install (editable)
pip install -e .

# Run all tests
pytest tests/

# Run a single test file
pytest tests/test_analyzer.py

# Run a single test by name
pytest tests/test_analyzer.py::test_function_name

# CLI entry points (after install)
init_app --confdir <config> --logdir <log> --outdir <output> --dbdir <db>
install_db --confdir <config>
check_install --confdir <config>
```

There is no Makefile, linter configuration, or CI pipeline — development relies on pytest directly.

## Architecture

The system uses a **manual dependency injection** pattern. `VEBenchmarkContainer` (`aigct/container.py`) is the single wiring point — it reads `config/aigct.yaml`, constructs all repositories and services, and exposes them as properties. All production code and tests instantiate one container and access components through it.

### Component Layers

```
VEBenchmarkContainer
  ├── query_mgr   → VEBenchmarkQueryMgr    # query variants, tasks, VEPs, filters
  ├── analyzer    → VEAnalyzer             # compute ROC/PR/MWU/calibration metrics
  ├── reporter    → VEAnalysisReporter     # format text summaries
  ├── plotter     → VEAnalysisPlotter      # matplotlib charts (ROC, PR, MWU, calibration)
  ├── exporter    → VEAnalysisExporter     # write results to timestamped CSV dirs
  └── data_validator → VEDataValidator     # check PK uniqueness across all tables
```

Repositories sit below the services and are accessed only through the container's private attributes (`_label_repo`, `_score_repo`, etc.). They are backed by CSV flat files, not a database.

### Data Access

`RepoSessionContext` (`aigct/repository.py`) manages table definitions and resolves file paths. It is shared across all six repositories.

Two singleton cache classes (`DataCache`, `VariantEffectLabelCache`) — both extending `ParameterizedSingleton` — hold lazy-loaded DataFrames in memory. The cache key includes task/VEP parameters, so the same CSV is only read once per unique parameter combination.

The flat-file schema has a global table set in `data/` plus per-task subdirectories (`data/CANCER/`, `data/ADRD/`, etc.) containing `variant_effect_label.csv` and `variant_effect_score.csv`.

### Domain Model

Key dataclasses in `aigct/model.py`:
- `VEQueryCriteria` — encapsulates all filter parameters (genes, variant IDs, allele frequency range, named filters)
- `VEAnalysisResult` — holds computed metrics (ROC-AUC, PR-AUC, MWU) per VEP
- `VEAnalysisCalibrationResult` — holds calibration curve data and binned scores
- `VariantId` — typed (assembly, chrom, pos, ref, alt) primary key for variants

### Configuration

`aigct/config/aigct.yaml.sample` is the template. The live config file is `config/aigct.yaml` (not committed). `VEBenchmarkContainer` defaults to `./config/aigct.yaml` relative to the working directory. `aigct/util.py` contains the `Config` class that converts the YAML dict into a recursive attribute namespace.

## Testing

Tests live in `tests/`. `tests/conftest.py` instantiates a real `VEBenchmarkContainer` (no mocking) and provides fixtures for all components plus synthetic user scores derived by sampling and perturbing existing repository data.

The `task_vep_code` fixture is parametrized over all task/VEP combinations in `TEST_TASK_VEP_CODES` — tests using it run once per combination.

Tests require the live database (CSV files) to be present at the path configured in `config/aigct.yaml`.
