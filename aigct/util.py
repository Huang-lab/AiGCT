"""
Utility classes
"""

import threading


class ParameterizedSingleton:
    """
    Classes that wish to behave as threadsafe singletons can inherit from
    this class. To be used only by classes that have an initialization
    method that takes parameters. The class must implement an _init_once
    method instead of the normal __init__ method for initialization.
    It takes same parameters as __init__ method. By inheriting from this
    class all instantiations of the subclass will return the same instance.
    """
    _instance = None
    _lock = threading.Lock()

    def __new__(cls, *args, **kwargs):
        if cls._instance is None:
            with cls._lock:
                if cls._instance is None:
                    cls._instance = super(ParameterizedSingleton,
                                          cls).__new__(cls)
                    cls._instance._init_once(*args, **kwargs)
        return cls._instance


class Config:
    def __init__(self, config):
        for key, value in config.items():
            if isinstance(value, dict):
                value = Config(value)
            setattr(self, key, value)


def str_or_list_to_list(str_or_list: str | list[str]) -> list[str]:
    if type(str_or_list) is str:
        return [str_or_list]
    else:
        return str_or_list




