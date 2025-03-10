"""
Dependency injector utilities. Bogus implementation to mimick a proper 
DI provider. This class could be reimplemented if we decide to use
a proper one in the future.
"""

import os
from .date_util import now_str_compact


def new_line(out, num_nls: int = 1):
    nls = "\n" * num_nls
    out.write(nls)


def unique_file_name(*args, suffix=None):
    file_name = now_str_compact(os.path.join(*args))
    return file_name + "." + suffix if suffix is not None else file_name


def create_folder(folder):
    if not os.path.exists(folder):
        os.mkdir(folder)
