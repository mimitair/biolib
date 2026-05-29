import pandas as pd

def safeCastColumns(df: pd.DataFrame, mapping: dict) -> pd.DataFrame:
    """ Allows typecasting of columns in a dataframe without raising errors if the column does not exist, or type is incompatible.
        Input:
           - df: pd.DataFrame:
           - mapping: dict

         Returns:
           - pd.DataFrame:

            """
    for column_name, dtype in mapping.items():
        try:
            df[column_name] = df[column_name].astype(dtype)
        except KeyError as e:
            print(f"'{column_name}' could not be found. Skipping typecast.")
        except ValueError as e:
            print(f"Error converting '{column_name}' to '{dtype}'")

    return df
