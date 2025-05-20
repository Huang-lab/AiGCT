"""
Class to simulate a proper Dependency Injection container.
It could be reimplemented in the future if we decide to use
a proper one. The interface, however, would remain the same.
"""

from ..util import Config
from ..repository import (
    RepoSessionContext,
    TABLE_DEFS
)
from .repo_loader import RepositoryLoader

import yaml
import os


class VEETLContainer:

    def __init__(self, config_file: str = "./config/config.yaml"):
        with (open(config_file, "r") as
              conf_file):
            self.config = Config(yaml.safe_load(conf_file))
        self._repo_session_context = RepoSessionContext(
            self.config.repository.root_dir, TABLE_DEFS)
        self._loader = RepositoryLoader(self.config,
                                        self._repo_session_context)

    @property
    def loader(self):
        return self._loader
