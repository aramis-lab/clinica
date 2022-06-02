
import pandas as pd
from string import Template


MISSING_TERM_ERROR_MSG = Template(
    "Term ${term} from the design matrix is not in the columns of the "
    "provided TSV file. Please make sure that there is no typo."
)


def _is_categorical(df: pd.DataFrame, column: str) -> bool:
    """This function returns whether the provided DataFrame's column
    is categorical or not.

    .. note::
        There might be more clever ways to do that.

    Parameters
    ----------
    df : DataFrame to analyze.
    column : Name of the column to check.

    Returns
    -------
    True if the column contains categorical values, False otherwise.
    """
    if column not in df.columns:
        raise ValueError(MISSING_TERM_ERROR_MSG.safe_substitute(term=column))
    return not df[column].dtype.name.startswith("float")
