from string import Template
from typing import Optional

import pandas as pd
from brainstat.stats.SLM import SLM
from brainstat.stats.terms import FixedEffect

__all__ = ["build_model", "print_clusters", "check_column_in_df", "is_categorical"]


MISSING_TERM_ERROR_MSG = Template(
    "Term ${term} from the design matrix is not in the columns of the "
    "provided TSV file. Please make sure that there is no typo."
)


def print_clusters(model: SLM, threshold: float) -> None:
    """This function prints the results related to total number
    of clusters, as well as the significative clusters.

    Parameters
    ----------
    model : brainstat.stats.SLM
        Fitted SLM model.

    threshold : float
        Cluster defining threshold.
    """
    from clinica.utils.stream import cprint

    cprint("#" * 40)
    cprint("After correction (Cluster-wise Correction for Multiple Comparisons): ")
    df = model.P["clus"][0]
    cprint(df)
    cprint(f"Clusters found: {len(df)}")
    cprint(
        f"Significative clusters (after correction): {len(df[df['P'] <= threshold])}"
    )


def build_model(design_matrix: str, df: pd.DataFrame) -> FixedEffect:
    """Build a brainstat model from the design matrix in string format.

    This function assumes that the design matrix is formatted
    in the following way:

        1 + factor_1 + factor_2 + ...

    Or:

        factor_1 + factor_2 + ...

    in the latter case the intercept will be added automatically.

    Parameters
    ----------
    design_matrix : str
        Design matrix specified as a string.

    df : pd.DataFrame
        Subjects DataFrame.

    Returns
    -------
    model : FixedEffect
        BrainStats model.
    """
    from functools import reduce

    if len(design_matrix) == 0:
        raise ValueError("Design matrix cannot be empty.")
    if "+" in design_matrix:
        terms = [_.strip() for _ in design_matrix.split("+")]
    else:
        terms = [design_matrix.strip()]
    model = []
    for term in terms:
        # Intercept is automatically included in brainstat
        if term == "1":
            continue
        # Handles the interaction effects
        if "*" in term:
            sub_terms = [_.strip() for _ in term.split("*")]
            model_term = reduce(
                lambda x, y: x * y, [_build_model_term(_, df) for _ in sub_terms]
            )
        else:
            model_term = _build_model_term(term, df)
        model.append(model_term)
    return reduce(lambda x, y: x + y, model)


def _build_model_term(
    term: str,
    df: pd.DataFrame,
    add_intercept: Optional[bool] = True,
) -> FixedEffect:
    """Builds a BrainStats model term from the subjects
    DataFrame and a column name.

    Parameters
    ----------
    term : str
        The name of the column of the DataFrame to be used.

    df : pd.DataFrame
        The subjects DataFrame.

    add_intercept : bool
        If `True`, adds an intercept term.

    Returns
    -------
    FixedEffect :
        BrainStats model term.
    """
    return FixedEffect(df[term], add_intercept=add_intercept)


def check_column_in_df(df: pd.DataFrame, column: str) -> None:
    """Checks if the provided column name is in the provided DataFrame.
    Raises a ValueError if not.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to analyze.

    column : str
        Name of the column to check.
    """
    if column not in df.columns:
        raise ValueError(MISSING_TERM_ERROR_MSG.safe_substitute(term=column))


def is_categorical(df: pd.DataFrame, column: str) -> bool:
    """Returns `True` if the column is categorical and `False` otherwise.

    Parameters
    ----------
    df : pd.DataFrame
        The DataFrame to analyze.

    column : str
        The name of the column to check.

    Returns
    -------
    bool :
        `True` if the column contains categorical values, `False` otherwise.
    """
    return not df[column].dtype.name.startswith("float")
