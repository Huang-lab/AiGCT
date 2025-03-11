# AI Genomics CollecTive (AIGCT)

## Overview

AI Genomics CollecTive (AIGCT) is a platform for systematically 
evaluating ML/AI models of variant effects across the spectrum of 
genomics-based precision medicine.

## Development/Testing Installation Instructions

Use these instructions for installation of the development version for internal testing
rather than the installation instructions in the readthedocs documentation below.
For usage refer to the readthedocs documentation.

Download the python package distribution file from sharepoint:

https://mtsinai-my.sharepoint.com/:u:/r/personal/kuan-lin_huang_mssm_edu/Documents/Huang_lab/manuscripts/AIPrecisionGenomics/demo/aigct-0.1a1.dev3.tar.gz?csf=1&web=1&e=cgC93V

Assume <apptarfile> is the full path to the downloaded file.

Install tarfile.

    pip install <apptarfile> --upgrade

Create a config, log, output, and db directory. Let’s refer to them as <config>, <log>, <output>, <db> for now.

Run the following to initialize the app and download the database files.

    init_app --confdir <config> --logdir <log> --outdir <output> --dbdir <dbdir>

where: <config>, <log>, <output>, and <dbdir> are directories where to store config file, log files, analysis output files, and database files, respectively.

Download the database file from sharepoint:

https://mtsinai-my.sharepoint.com/:u:/r/personal/kuan-lin_huang_mssm_edu/Documents/Huang_lab/manuscripts/AIPrecisionGenomics/demo/repo_v0_1a1_dev3.tar.gz?csf=1&web=1&e=Ld5yg0

Assume <dbtarfile> is full path to downloaded file.

Extract contents:

    tar -xf <dbtarfile> -C  <dbdir>

where <dbdir> is the database directory specified for init_app above.


## Documentation Including Installation and Usage Instructions

**[readthedocs](https://aigct-dev.readthedocs.io/en/latest/index.html#)**
