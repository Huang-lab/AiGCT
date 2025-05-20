from datetime import datetime


def now_str_compact(prefix: str = ""):
    now = datetime.now()
    return now.strftime(prefix + "%Y%m%d%H%M%S")


def now_str_basic_format():
    now = datetime.now()
    return now.strftime("%Y/%m/%d %H:%M:%S")


