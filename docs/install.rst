Installation
============


Install python package::

    pip install aigct

Initialize the application::

    init_app –-confdir <config> --logdir <log> --outdir <output> --dbdir <dbdir>

Where <config>, <log>, <output>, <dbdir> are directories where to store
config file, log files, output files containing the results of a benchmarking
analysis, and database files, respectively.

Edit the config file, <config>/aigct.yaml, setting the source_url parameter to
the url of the database files. For the latest release it is:

``https://zenodo.org/records/19928885/files/aigct_repo.tar.gz?download=1``

Download and install the database files::

    install_db --confdir <config>

You should see .csv files in the <dbdir> directory.


