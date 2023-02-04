from os import PathLike
from typing import List, Optional


class MissingModsTracker:
    """Class used for tracking the number of missing modalities in a database.

    Attributes
    ----------
    missing : dict
        Dictionary of missing modalities.

    ses : list of str
        Sessions for which the tracker has information.
    """

    _default_modalities_tracked = ["dwi", "func", "fieldmap", "flair", "t1w"]

    def __init__(self, sessions: List[str], mod_list: Optional[List[str]] = None):
        modalities = mod_list or self._default_modalities_tracked
        if "session" not in modalities:
            modalities.insert(0, "session")
        self.missing = {s: {mod: 0 for mod in modalities} for s in sessions}

    @property
    def ses(self):
        return list(self.missing.keys())

    def add_missing_mod(self, session: str, modality: str) -> None:
        """Increase the number of missing files for the input modality.

        Parameters
        ----------
        session : str
            Name of the session for which to add a modality.

        modality : str
            The missing modality to add.
        """
        if session not in self.missing:
            raise ValueError(
                f"Session {session} was not provided to the MissingModsTracker constructor."
            )
        if modality not in self.missing[session]:
            raise ValueError(
                f"Modality {modality} is not tracked by this instance of MissingModsTracker."
            )
        self.missing[session][modality] += 1

    def increase_missing_ses(self, session: str) -> None:
        """Increase the session number.

        Parameters
        ----------
        session : str
            Name of the session.
        """
        if session not in self.missing:
            raise ValueError(
                f"Session {session} was not provided to the MissingModsTracker constructor."
            )
        self.missing[session]["session"] += 1

    def get_missing_list(self) -> dict:
        """Return the dictionary of missing modalities.

        Returns
        -------
        dict :
             The dictionary containing the list of missing files.
        """
        return self.missing


def sort_session_list(session_list: List[str]) -> List[str]:
    """Sorts the list of session IDs provided based on their session number.

    Parameters
    ------------
    session_list : list[str]
        List of session IDs to sort.

    Returns
    --------
    list[str] :
        Sorted list of session IDs.

    Examples
    --------
    >>>sort_session_list(["ses-M000", "ses-M006", "ses-M012", "ses-M024", "ses-M048", "ses-M003"])
    ["ses-M000", "ses-M003", "ses-M006", "ses-M012", "ses-M024", "ses-M048"]
    >>>sort_session_list(["ses-M0", "ses-M6", "ses-M12", "ses-M24", "ses-M48", "ses-M3"])
    ["ses-M0", "ses-M3", "ses-M6", "ses-M12", "ses-M24", "ses-M48"]
    """
    prefix_length = len("ses-M")
    session_idx = [
        ((session[prefix_length:]), int(session[prefix_length:]))
        for session in session_list
    ]
    session_idx.sort(key=lambda x: x[1])
    return [f"ses-M{session[0]}" for session in session_idx]


def write_statistics(
    summary_file: PathLike,
    num_subjs: int,
    ses_avail: List[str],
    mmt: MissingModsTracker,
) -> None:
    """Write statistics file.

    Write to a given input file statistics about missing files and modalities in a dataset.
    This method takes in input a MissingModsTracker object (mmt)
    that contains the number of missing modalities for each session and the
    number of missing sessions for each subject.

    Parameters
    ----------
    summary_file : PathLike
        Path of the output file where statistics should be written.

    num_subjs : int
        Number of subjects.

    ses_avail : list of str
        List of sessions available.

    mmt : MissingModsTracker
        Instance of MissingModsTracker.
    """
    with open(summary_file, "w") as fp:
        fp.write(compute_statistics(num_subjs, ses_avail, mmt))


def compute_statistics(
    num_subjs: int, ses_avail: List[str], mmt: MissingModsTracker
) -> str:
    """Compute statistics and return them as a string for printing/writing.

    Statistics are about missing files and modalities in a dataset.
    This method takes in input a MissingModsTracker object (mmt)
    that contains the number of missing modalities for each session and the
    number of missing sessions for each subject.

    Parameters
    ----------
    num_subjs : int
        Number of subjects.

    ses_avail : list of str
        List of sessions available.

    mmt : MissingModsTracker
        Instance of MissingModsTracker.

    Returns
    -------
    summary : str
        The statistics formatted as a summary string.
    """
    missing_list = mmt.get_missing_list()
    ses_avail = sort_session_list(ses_avail)
    ses_founds = {ses: num_subjs - missing_list[ses]["session"] for ses in ses_avail}
    summary = "\n".join(
        [
            "*" * 46,
            f"Number of subjects converted: {num_subjs}",
            f"Sessions available: {ses_avail}\n",
            "\n".join(
                [
                    f"Number of sessions {ses} found: {ses_found} ({ses_found * 100 / num_subjs}%)\n"
                    for ses, ses_found in ses_founds.items()
                ]
            ),
            "*" * 46 + "\n\n" + "Number of missing modalities for each session:\n",
        ]
    )
    for ses in ses_avail:
        summary += "\n" + ses + "\n"
        for mod in missing_list[ses]:
            if mod != "session":
                num_miss_mod = missing_list[ses][mod]
                percentage_missing = round(
                    (
                        num_miss_mod
                        * 100
                        / float(num_subjs - missing_list[ses]["session"])
                    ),
                    2,
                )
                summary += f"{mod}: {num_miss_mod} ({percentage_missing}%) \n"

    return summary


def write_longitudinal_analysis(
    summary_file: PathLike,
    bids_dir: PathLike,
    out_dir: PathLike,
    ses_avail: List[str],
    out_file_name: str,
) -> None:
    """Write to a given input file statistics about the present modalities and diagnoses in a dataset for each session.

    Parameters
    ----------
    summary_file : PathLike
        Path to the file where statistics should be written.

    bids_dir : PathLike
        Path to the BIDS directory.

    out_dir : PathLike
        Path to the output directory of the check-missing-modality pipeline.

    ses_avail : list of str
        List of sessions available.

    out_file_name : str
        String that replaces the default prefix ('missing_mods_') in
        the name of all the output files.
    """
    with open(summary_file, "w") as fp:
        fp.write(
            compute_longitudinal_analysis(bids_dir, out_dir, ses_avail, out_file_name)
        )


def compute_longitudinal_analysis(
    bids_dir: PathLike,
    out_dir: PathLike,
    ses_avail: List[str],
    out_file_name: str,
) -> str:
    """Compute statistics about the present modalities and diagnoses in a dataset for each session.

    Parameters
    ----------
    bids_dir : PathLike
        Path to the BIDS directory.

    out_dir : PathLike
        Path to the output directory of the check-missing-modality pipeline.

    ses_avail : list of str
        List of sessions available.

    out_file_name : str
        String that replaces the default prefix ('missing_mods_') in
        the name of all the output files.

    Returns
    -------
    summary : str
        The summary analysis as a string.
    """
    from collections import Counter
    from pathlib import Path

    import pandas as pd

    out_dir = Path(out_dir)
    bids_dir = Path(bids_dir)
    ses_avail = sort_session_list(ses_avail)
    summary = "\n\n".join(
        ["*" * 46, f"Number of present diagnoses and modalities for each session:\n"]
    )
    for ses in ses_avail:
        ses_df = pd.read_csv(
            out_dir / (out_file_name + ses + ".tsv"), sep="\t"
        ).set_index("participant_id")
        mods_avail = ses_df.columns.values
        mod_dict = dict()
        for mod in mods_avail:
            diagnosis_counter = Counter()
            subjects_avail = ses_df[ses_df[mod] == 1].index.values
            for subject in subjects_avail:
                subj_tsv_path = bids_dir / subject / f"{subject}_sessions.tsv"
                subj_df = pd.read_csv(subj_tsv_path, sep="\t").set_index("session_id")
                diagnosis = "missing"
                if (
                    ses in subj_df.index.values
                    and "diagnosis" in subj_df.columns.values
                ):
                    diagnosis = subj_df.loc[ses, "diagnosis"]
                    diagnosis = diagnosis if isinstance(diagnosis, str) else "n/a"
                diagnosis_counter.update([diagnosis])
            mod_dict[mod] = dict(diagnosis_counter)
        summary += f"{ses}\n{compute_table(mod_dict)}\n\n"
    return summary


def compute_table(mod_dict: dict) -> str:
    """Builds a table, encoded as a string, describing the
    modalities in the given dictionary.

    Parameters
    ----------
    mod_dict : dict
        Dictionary of modalities for which to compute the table.

    Returns
    -------
    table : str
        The summary table as a string.
    """
    diagnoses = sorted(list(set().union(*mod_dict.values())))
    table = "\t" + "\t| ".join(diagnoses)
    table += "\n" + "-" * 8 * (len(diagnoses) + 2) + "\n"
    for mod_key, mod_value in mod_dict.items():
        table += f"{mod_key}"
        if len(mod_key) < 8:
            table += "\t"
        table += "\t| ".join([str(mod_value.get(d, "0")) for d in diagnoses]) + "\n"
    return table


def viscode_to_session(viscode: str) -> str:
    """Replace the session label 'bl' with 'M000' or capitalize the session name passed as input.

    Parameters
    ----------
    viscode: str
        The name of the session.

    Returns
    -------
    str:
        "M000" if the session is the baseline session. Otherwise returns the original session name capitalized.
    """
    if viscode in {"bl", "m0"}:
        return "ses-M000"
    else:
        return "ses-" + f"M{(int(viscode[1:])):03d}"
