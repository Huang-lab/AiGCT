import pandas as pd
from .util import str_or_list_to_list

RELATIONAL_OPERATORS = ["==", "!=", ">", "<", ">=", "<="]


def filter_dataframe_by_list(data_frame: pd.DataFrame,
                             filter_list: pd.DataFrame | list[str] | str |
                             pd.Series,
                             df_merge_columns: list[str] | str,
                             filter_col_name_map: dict = None,
                             in_list: bool = True
                             ) -> pd.DataFrame:
    if len(data_frame) == 0:
        return data_frame
    if len(filter_list) == 0:
        return data_frame[:0].copy()
    df_merge_columns = str_or_list_to_list(df_merge_columns)
    if filter_col_name_map is None:
        filter_merge_columns = df_merge_columns
    else:
        filter_merge_columns = [filter_col_name_map[merge_col] for
                                merge_col in df_merge_columns]
    if type(filter_list) is pd.DataFrame:
        filter_df = filter_list[filter_merge_columns].drop_duplicates()
    else:
        if len(df_merge_columns) > 1:
            raise Exception("Cannot filter a dataframe by more than " +
                            f"one column: {df_merge_columns}: " +
                            "when filtering by a list rather than by a " +
                            "dataframe.")
        if type(filter_list) is str:
            filter_list = [filter_list]
        filter_df = pd.DataFrame({df_merge_columns[0]: list(set(filter_list))})
    if in_list:
        return data_frame.merge(filter_df[filter_merge_columns],
                                left_on=df_merge_columns,
                                right_on=filter_merge_columns,
                                how="inner")[data_frame.columns]
    data_frame_ret = data_frame.merge(filter_df[filter_merge_columns],
                                      left_on=df_merge_columns,
                                      right_on=filter_merge_columns,
                                      how="left", indicator=True)
    return data_frame_ret.query("_merge == 'left_only'").drop(
        columns="_merge")[data_frame.columns]


def build_dataframe_where_clause(where_params: dict) -> str:
    clause: str = ""
    for column, operator_value in where_params.items():
        if operator_value[1] is None:
            continue
        operator = operator_value[0]
        if operator not in RELATIONAL_OPERATORS:
            raise Exception(
                f"{operator} is not a legal query relational operator. " +
                f"It must be one of {RELATIONAL_OPERATORS}")
        if clause != "":
            clause += " "
        clause = f"{clause}{column} {operator} "
        if type(operator_value[1]) is str:
            val = "\"" + operator_value[1] + "\""
        else:
            val = str(operator_value[1])
        clause += val
    return clause


def merge_by_chunks(left_df: pd.DataFrame, right_df: pd.DataFrame,
                    chunk_size: int = 500, **kwargs):
    chunks = []
    for i in range(0, len(left_df), chunk_size):
        chunk = left_df.iloc[i:i+chunk_size].merge(right_df, **kwargs)
        chunks.append(chunk)
    return pd.concat(chunks)
