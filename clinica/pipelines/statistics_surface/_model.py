from functools import reduce

import pandas as pd
from brainstat.stats.terms import FixedEffect

from ._utils import MISSING_TERM_ERROR_MSG, _is_categorical


def _build_model(design_matrix: str, df: pd.DataFrame):
    """Build a brainstat model from the design matrix in
    string format.

    This function assumes that the design matrix is formatted
    in the following way:

        1 + factor_1 + factor_2 + ...

    Or:

        factor_1 + factor_2 + ...

    in the latter case the intercept will be added automatically.

    Parameters
    ----------
    design_matrix : Design matrix specified as a string.
    df : Subjects DataFrame.

    Returns
    -------
    model : BrainStats model.
    """
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
    if len(model) == 1:
        return model[0]
    return reduce(lambda x, y: x + y, model)


def _build_model_term(term: str, df: pd.DataFrame) -> FixedEffect:
    """Builds a BrainStats model term from the subjects
    DataFrame and a column name.

    Parameters
    ----------
    term : Name of the column of the DataFrame to be used.
    df : Subjects DataFrame.

    Returns
    -------
    BrainStats FixedEffect.
    """
    if term not in df.columns:
        raise ValueError(MISSING_TERM_ERROR_MSG.safe_substitute(term=term))
    return FixedEffect(df[term])
