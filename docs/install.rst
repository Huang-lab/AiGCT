Installation
============


Install python package::

    pip install aigct

Initialize the application::

    init_app –-confdir <config> --logdir <log> --outdir <output> --dbdir <dbdir>

Where <config>, <log>, <output>, <dbdir> are directories where to store
config file, log files, output files containing the results of a benchmarking
analysis, and database files, respectively.

Download and install the database files::

    install_db --confdir <config>

You should see .csv files in the <dbdir> directory.


