"""Session-level usability filtering and cprint-based logging."""

import pandas as pd

from ._scan_classification import USABLE_ACTIONS

SESSION_KEY = ["Subject_ID", "Session_ID"]


def split_usable_sessions(
    df: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Split the scan inventory into usable and unusable sessions.

    A session is *usable* when at least one of its runs has an Action that
    belongs to ``USABLE_ACTIONS``.  Sessions where *no* run is usable are
    returned as the unusable set.

    Parameters
    ----------
    df:
        Full scan inventory DataFrame (one row per run).

    Returns
    -------
    (usable_df, unusable_df)
        Both DataFrames share the same columns as the input.
    """
    has_usable = df.groupby(SESSION_KEY)["Action"].transform(
        lambda actions: actions.isin(USABLE_ACTIONS).any()
    )
    return df[has_usable].copy(), df[~has_usable].copy()


def log_session_summary(usable_df: pd.DataFrame, unusable_df: pd.DataFrame) -> None:
    """Log a session-level summary using ``cprint``.

    Logs:
    - Total number of sessions evaluated
    - Number of usable sessions (at least one run fit for SUVR)
    - Number of unusable sessions
    - Subject_ID and Session_ID for each unusable session
    """
    from clinica.utils.stream import cprint

    total = (
        usable_df.groupby(SESSION_KEY).ngroups
        + unusable_df.groupby(SESSION_KEY).ngroups
    )
    n_usable = usable_df.groupby(SESSION_KEY).ngroups
    n_unusable = unusable_df.groupby(SESSION_KEY).ngroups

    cprint(f"Sessions evaluated:             {total}", lvl="info")
    cprint(f"Sessions with a usable scan:    {n_usable}", lvl="info")
    cprint(
        f"Sessions with NO usable scan:   {n_unusable}",
        lvl="warning" if n_unusable else "info",
    )

    if n_unusable:
        lines = ["Sessions with no usable AV1451 scan:"]
        for _, grp in unusable_df.groupby(SESSION_KEY):
            row = grp.iloc[0]
            lines.append(f"  {row['Subject_ID']}  {row['Session_ID']}")
        cprint("\n".join(lines), lvl="warning")
