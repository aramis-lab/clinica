def sort_session_list(session_list):
    session_idx = [int(session[5:]) for session in session_list]
    session_idx.sort()
    session_id_list = []
    for session in session_idx:
        if session < 10:
            session_id_list.append(f"ses-M0{session}")
        else:
            session_id_list.append(f"ses-M{session}")

    return session_id_list


def print_statistics(summary_file, num_subjs, ses_aval, mmt):
    """Write statistics file.

    Print to a given input file statistics about missing files and modalities
    in a dataset.  This metod takes in input a MissingModsTracker object (mmt)
    that contains the number of missing modalities for each session and the
    number of missing sessions for each subject.

    Args:
        summary_file: path of the output file where write.
        num_subjs: number of subjects.
        ses_aval: list of sessions available.
        mmt: object MissingModsTracker
    """
    missing_list = mmt.get_missing_list()
    ses_aval = sort_session_list(ses_aval)
    summary_file.write("**********************************************\n")
    summary_file.write(f"Number of subjects converted: {num_subjs}\n")
    summary_file.write(f"Sessions available: {ses_aval}\n")

    ses_aval = sort_session_list(ses_aval)

    for ses in ses_aval:
        ses_miss = missing_list[ses]["session"]
        ses_found = num_subjs - ses_miss
        perc_ses_found = ses_found * 100 / num_subjs
        summary_file.write(
            f"Number of sessions {ses} found: {ses_found} ({perc_ses_found}%)\n"
        )

    summary_file.write("**********************************************\n\n")
    summary_file.write("Number of missing modalities for each session:\n")

    for ses in ses_aval:
        summary_file.write("\n" + ses + "\n")
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
                summary_file.write(f"{mod}: {num_miss_mod} ({percentage_missing}%) \n")


def print_longitudinal_analysis(
    summary_file, bids_dir, out_dir, ses_aval, out_file_name
):
    """Print to a given input file statistics about the present modalities and diagnoses in a dataset for each session.

    Args:
        summary_file: path of the output file where write.
        bids_dir: path to the BIDS directory
        out_dir: output_dir of the check-missing-modality pipeline.
        ses_aval: list of sessions available.
        out_file_name: string that replace the default prefix ('missing_mods_') in the name of all the output files
    """
    from os import path

    import pandas as pd

    ses_aval = sort_session_list(ses_aval)

    summary_file.write("**********************************************\n\n")
    summary_file.write("Number of present diagnoses and modalities for each session:\n")

    for ses in ses_aval:
        ses_df = pd.read_csv(
            path.join(out_dir, out_file_name + ses + ".tsv"), sep="\t"
        ).set_index("participant_id")
        mods_avail = ses_df.columns.values
        mod_dict = dict()
        for mod in mods_avail:
            diagnosis_dict = dict()
            subjects_avail = ses_df[ses_df[mod] == 1].index.values
            for subject in subjects_avail:
                subj_tsv_path = path.join(bids_dir, subject, f"{subject}_sessions.tsv")
                subj_df = pd.read_csv(subj_tsv_path, sep="\t").set_index("session_id")
                if (
                    ses in subj_df.index.values
                    and "diagnosis" in subj_df.columns.values
                ):
                    diagnosis = subj_df.loc[ses, "diagnosis"]
                    if isinstance(diagnosis, str):
                        increment_dict(diagnosis_dict, diagnosis)
                    else:
                        increment_dict(diagnosis_dict, "n/a")
                else:
                    increment_dict(diagnosis_dict, "missing")

            mod_dict[mod] = diagnosis_dict

        summary_file.write(f"{ses}\n")
        print_table(summary_file, mod_dict)
        summary_file.write("\n\n")


def increment_dict(dictionnary, key):
    if key not in dictionnary:
        dictionnary[key] = 1
    else:
        dictionnary[key] += 1


def print_table(summary_file, double_dict):

    # Find all keys at the second level of the dictionnary
    diagnoses = set()
    mods = double_dict.keys()
    for mod in mods:
        diagnoses = diagnoses | set(double_dict[mod].keys())

    diagnoses = list(diagnoses)
    diagnoses.sort()

    summary_file.write("\t")
    for diag in diagnoses:
        summary_file.write(f"\t| {diag}")
    summary_file.write("\n" + "-" * 8 * (len(diagnoses) + 2) + "\n")

    for mod in double_dict.keys():
        summary_file.write(f"{mod}")
        if len(mod) < 8:
            summary_file.write("\t")
        for diag in diagnoses:
            if diag not in double_dict[mod]:
                summary_file.write("\t| 0")
            else:
                summary_file.write(f"\t| {double_dict[mod][diag]}")
        summary_file.write("\n")


def has_one_index(index_list):
    if len(index_list) == 1:
        return index_list[0]
    if len(index_list) == 0:
        return -1
    if len(index_list) > 1:
        raise ValueError("Multiple indexes found")


class MissingModsTracker:
    """Class used for tracking the number of missing modalities in a database."""

    def __init__(self, ses, mod_list=None):
        self.missing = {}
        self.ses = ses
        if mod_list:
            for s in ses:
                self.missing.update({s: {"session": 0}})
                for mod in mod_list:
                    self.missing[s].update({mod: 0})
        else:
            for s in ses:
                self.missing.update(
                    {
                        s: {
                            "session": 0,
                            "dwi": 0,
                            "func": 0,
                            "fieldmap": 0,
                            "flair": 0,
                            "t1w": 0,
                        }
                    }
                )

    def add_missing_mod(self, ses, mod):
        """Increase the number of missing files for the input modality.

        Args:
            ses: name of the session
            mod: modality missing
        """
        self.missing[ses][mod] += 1

    def increase_missing_ses(self, ses):
        self.missing[ses]["session"] += 1

    def get_missing_list(self):
        """Return the hash map mods_missing.

        Returns:
             The hash map containing the list of missing files.
        """
        return self.missing
