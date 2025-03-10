from pandas.api.types import is_string_dtype
from .model import VEQueryCriteria


def validate_query_criteria(qry: VEQueryCriteria):
    if (qry is not None and qry.variant_ids is not None and
            len(qry.variant_ids) > 0):
        if not is_string_dtype(qry.variant_ids['CHROMOSOME']):
            raise Exception("Values in CHROMOSOME column must all be strings")
