# coding: utf8

"""Data handling scripts."""

import click


def compute_default_filename(out_path):
    from os import path

    abspath = path.abspath(out_path)
    # If given path is a directory, append a filename
    if path.isdir(abspath):
        tsv_path = path.join(out_path, "merge.tsv")
    elif "." not in path.basename(abspath):
        tsv_path = f"{out_path}.tsv"
    else:
        if path.splitext(out_path)[1] != ".tsv":
            raise TypeError("Output path extension must be tsv.")
        tsv_path = out_path

    return tsv_path


def create_merge_file(
    bids_dir,
    out_tsv,
    caps_dir=None,
    tsv_file=None,
    pipelines=None,
    ignore_scan_files=False,
    ignore_sessions_files=False,
    **kwargs,
):
    """Merge all the TSV files containing clinical data of a BIDS compliant dataset and store the result inside a TSV file.

    Args:
        bids_dir: path to the BIDS folder
        out_tsv: path to the output tsv file
        caps_dir: path to the CAPS folder (optional)
        tsv_file: TSV file containing the subjects with their sessions (optional)
        ignore_scan_files: If True the information related to scans is not read (optional)
        ignore_sessions_files: If True the information related to sessions and scans is not read (optional)
        pipelines: when adding CAPS information, indicates the pipelines that will be merged (optional)
    """
    import json
    import os
    from os import path

    import numpy as np
    import pandas as pd

    from clinica.utils.participant import get_subject_session_list
    from clinica.utils.stream import cprint

    from .pipeline_handling import DatasetError

    if caps_dir is not None:
        if not path.isdir(caps_dir):
            raise IOError("The path to the CAPS directory is wrong")

    if not os.path.isfile(path.join(bids_dir, "participants.tsv")):
        raise IOError("participants.tsv not found in the specified BIDS directory")
    participants_df = pd.read_csv(path.join(bids_dir, "participants.tsv"), sep="\t")

    sessions, subjects = get_subject_session_list(
        bids_dir, ss_file=tsv_file, use_session_tsv=True
    )
    sub_ses_df = pd.DataFrame(
        [[subject, session] for subject, session in zip(subjects, sessions)],
        columns=["participant_id", "session_id"],
    )
    sub_ses_df.set_index(["participant_id", "session_id"], inplace=True)

    out_path = compute_default_filename(out_tsv)
    out_dir = path.dirname(out_path)
    if len(out_dir) > 0:
        os.makedirs(out_dir, exist_ok=True)

    merged_df = pd.DataFrame(columns=participants_df.columns.values)

    # BIDS part
    for subject, subject_df in sub_ses_df.groupby(level=0):
        sub_path = path.join(bids_dir, subject)
        row_participant_df = participants_df[
            participants_df["participant_id"] == subject
        ]
        row_participant_df.reset_index(inplace=True, drop=True)
        if len(row_participant_df) == 0:
            cprint(
                msg=f"Participant {subject} does not exist in participants.tsv",
                lvl="warning",
            )
            row_participant_df = pd.DataFrame([[subject]], columns=["participant_id"])

        if ignore_sessions_files:
            for _, session in subject_df.index.values:
                row_session_df = pd.DataFrame([[session]], columns=["session_id"])

                row_df = pd.concat([row_participant_df, row_session_df], axis=1)
                merged_df = merged_df.append(row_df)

        else:
            sessions_df = pd.read_csv(
                path.join(sub_path, f"{subject}_sessions.tsv"), sep="\t"
            )

            for _, session in subject_df.index.values:
                row_session_df = sessions_df[sessions_df.session_id == session]
                row_session_df.reset_index(inplace=True, drop=True)
                if len(row_session_df) == 0:
                    raise DatasetError(
                        sessions_df.loc[0, "session_id"] + " / " + session
                    )

                # Read scans TSV files
                scan_path = path.join(
                    bids_dir,
                    subject,
                    session,
                    f"{subject}_{session}_scans.tsv",
                )
                if path.isfile(scan_path) and not ignore_scan_files:
                    scans_dict = dict()
                    scans_df = pd.read_csv(scan_path, sep="\t")
                    for idx in scans_df.index.values:
                        filepath = scans_df.loc[idx, "filename"]
                        if filepath.endswith(".nii.gz"):
                            filename = path.basename(filepath).split(".")[0]
                            modality = "_".join(filename.split("_")[2::])
                            for col in scans_df.columns.values:
                                if col == "filename":
                                    pass
                                else:
                                    value = scans_df.loc[idx, col]
                                    new_col_name = f"{modality}_{col}"
                                    scans_dict.update({new_col_name: value})
                            json_path = path.join(
                                bids_dir,
                                subject,
                                session,
                                filepath.split(".")[0] + ".json",
                            )
                            if path.exists(json_path):
                                with open(json_path, "r") as f:
                                    json_dict = json.load(f)
                                for key, value in json_dict.items():
                                    new_col_name = f"{modality}_{key}"
                                    scans_dict.update({new_col_name: value})
                    scans_dict = {
                        str(key): str(value) for key, value in scans_dict.items()
                    }
                    row_scans_df = pd.DataFrame(scans_dict, index=[0])
                else:
                    row_scans_df = pd.DataFrame()

                row_df = pd.concat(
                    [row_participant_df, row_session_df, row_scans_df], axis=1
                )
                merged_df = merged_df.append(row_df)

    # Put participant_id and session_id first
    col_list = merged_df.columns.values.tolist()
    col_list.insert(0, col_list.pop(col_list.index("participant_id")))
    col_list.insert(1, col_list.pop(col_list.index("session_id")))
    merged_df = merged_df[col_list]

    tmp = merged_df.select_dtypes(include=[np.number])
    # Round numeric values in dataframe to 6 decimal values
    merged_df.loc[:, tmp.columns] = np.round(tmp, 6)
    merged_df.to_csv(out_path, sep="\t", index=False)
    cprint("End of BIDS information merge.", lvl="debug")

    merged_df.reset_index(drop=True, inplace=True)

    # CAPS
    if caps_dir is not None:
        # Call the different pipelines
        from .pipeline_handling import (
            pet_volume_pipeline,
            t1_freesurfer_pipeline,
            t1_volume_pipeline,
        )

        pipeline_options = {
            "t1-volume": t1_volume_pipeline,
            "pet-volume": pet_volume_pipeline,
            "t1-freesurfer": t1_freesurfer_pipeline,
        }
        merged_summary_df = pd.DataFrame()
        if not pipelines:
            for pipeline_name, pipeline_fn in pipeline_options.items():
                merged_df, summary_df = pipeline_fn(caps_dir, merged_df, **kwargs)
                if summary_df is not None and not summary_df.empty:
                    merged_summary_df = pd.concat([merged_summary_df, summary_df])

                if summary_df is None or summary_df.empty:
                    cprint(
                        f"{pipeline_name} outputs were not found in the CAPS folder."
                    )
        else:
            for pipeline in pipelines:
                merged_df, summary_df = pipeline_options[pipeline](
                    caps_dir, merged_df, **kwargs
                )
                merged_summary_df = pd.concat([merged_summary_df, summary_df])

        n_atlas = len(merged_summary_df)
        if n_atlas == 0:
            raise FileNotFoundError(
                "No outputs were found for any pipeline in the CAPS folder. "
                "The output only contains BIDS information."
            )
        columns = merged_df.columns.values.tolist()
        merged_summary_df.reset_index(inplace=True, drop=True)
        for idx in merged_summary_df.index:
            first_column_name = merged_summary_df.loc[idx, "first_column_name"]
            last_column_name = merged_summary_df.loc[idx, "last_column_name"]
            merged_summary_df.loc[idx, "first_column_index"] = columns.index(
                first_column_name
            )
            merged_summary_df.loc[idx, "last_column_index"] = columns.index(
                last_column_name
            )

        summary_path = path.splitext(out_path)[0] + "_summary.tsv"
        merged_summary_df.to_csv(summary_path, sep="\t", index=False)

        tmp = merged_df.select_dtypes(include=[np.number])
        # Round numeric values in dataframe to 12 floating point values
        merged_df.loc[:, tmp.columns] = np.round(tmp, 12)
        merged_df.to_csv(out_path, sep="\t")
        cprint("End of CAPS information merge.", lvl="debug")


def find_mods_and_sess(bids_dir):
    """Find all the modalities and sessions available for a given BIDS dataset.

    Args:
        bids_dir: path to the BIDS dataset

    Returns:
        mods_dict: a dictionary that stores the sessions and modalities found and has the following structure.
    Example:
    {
        'sessions': ['ses-M00', 'ses-M18'],
        'fmap': ['fmap'],
        'anat': ['flair', 't1w'],
        'func': ['func_task-rest'],
        'dwi': ['dwi']
    }

    """
    import os
    from glob import glob
    from os import path

    mods_dict = {}
    mods_list = []
    subjects_paths_lists = glob(path.join(bids_dir, "*sub-*"))

    for sub_path in subjects_paths_lists:
        ses_paths = glob(path.join(sub_path, "*ses-*"))
        for session in ses_paths:
            ses_name = session.split(os.sep)[-1]
            mods_avail = []
            if "sessions" in mods_dict:
                if ses_name not in mods_dict["sessions"]:
                    mods_dict["sessions"].append(ses_name)
            else:
                mods_dict.update({"sessions": [ses_name]})
            mods_paths_folders = glob(path.join(session, "*/"))

            for p in mods_paths_folders:
                p = p[:-1]
                mods_avail.append(p.split("/").pop())
            if "func" in mods_avail:
                list_funcs_paths = glob(path.join(session, "func", "*bold.nii.gz"))
                for func_path in list_funcs_paths:
                    func_name = func_path.split(os.sep)[-1]
                    func_name_tokens = func_name.split("_")
                    func_task = func_name_tokens[2]
                    if "func" in mods_dict:
                        if "func_" + func_task not in mods_dict["func"]:
                            mods_dict["func"].append("func_" + func_task)
                    else:
                        mods_dict.update({"func": ["func_" + func_task]})

                    if "func_" + func_task not in mods_list:
                        mods_list.append("func_" + func_task)

            if "dwi" in mods_avail:
                if "dwi" not in mods_dict:
                    mods_dict.update({"dwi": ["dwi"]})
                if "dwi" not in mods_list:
                    mods_list.append("dwi")

            if "fmap" in mods_avail:
                if "fmap" not in mods_dict:
                    mods_dict.update({"fmap": ["fmap"]})
                if "fmap" not in mods_list:
                    mods_list.append("fmap")

            if "pet" in mods_avail:
                list_pet_paths = glob(path.join(session, "pet", "*pet.nii.gz"))
                for pet_path in list_pet_paths:
                    pet_name = pet_path.split(os.sep)[-1].split(".")[0]
                    pet_name_tokens = pet_name.split("_")
                    pet_acq = pet_name_tokens[3]
                    if "pet" in mods_dict:
                        if "pet_" + pet_acq not in mods_dict["pet"]:
                            mods_dict["pet"].append("pet_" + pet_acq)
                    else:
                        mods_dict.update({"pet": ["pet_" + pet_acq]})

                    if "pet_" + pet_acq not in mods_list:
                        mods_list.append("pet_" + pet_acq)

            if "anat" in mods_avail:
                anat_files_paths = glob(path.join(session, "anat", "*"))

                for anat_file in anat_files_paths:
                    anat_name = anat_file.split(os.sep)[-1]

                    # Extract the name of the file without the extension
                    if ".nii.gz" in anat_name:
                        anat_name = anat_name.replace(".nii.gz", "")
                        anat_ext = "nii.gz"
                    else:
                        anat_name = os.path.splitext(anat_name.split(os.sep)[-1])[0]
                        anat_ext = os.path.splitext(anat_name.split(os.sep)[-1])[1]

                    if anat_ext != "json":
                        file_parts = anat_name.split("_")
                        anat_type = str.lower(file_parts[len(file_parts) - 1])
                        if "anat" in mods_dict:
                            if anat_type not in mods_dict["anat"]:
                                anat_aval = mods_dict["anat"]
                                anat_aval.append(anat_type)
                                mods_dict.update({"anat": anat_aval})
                        else:
                            mods_dict.update({"anat": [anat_type]})

                        if anat_type not in mods_list:
                            mods_list.append(anat_type)

    return mods_dict


def compute_missing_processing(bids_dir, caps_dir, out_file):
    """
    Compute the list of missing processing for each subject in a CAPS compliant dataset

    Args:
        bids_dir: path to the BIDS directory.
        caps_dir: path to the CAPS directory.
        out_file: path to the output file (filename included).
    """
    from glob import glob
    from os import listdir, path, sep

    import pandas as pd

    if path.exists(path.join(caps_dir, "groups")):
        groups = listdir(path.join(caps_dir, "groups"))
    else:
        groups = list()
    output_df = pd.DataFrame()

    # Retrieve pet tracers avail
    mods_and_sess = find_mods_and_sess(bids_dir)
    mods_and_sess.pop("sessions")
    mods_avail_dict = mods_and_sess
    trc_avail = [
        j.split("_")[1] for i in mods_avail_dict.values() for j in i if "pet" in j
    ]
    print(trc_avail)

    subjects_paths = glob(path.join(caps_dir, "subjects", "sub-*"))
    for subject_path in subjects_paths:
        participant_id = subject_path.split(sep)[-1]
        sessions_paths = glob(path.join(subject_path, "ses-*"))
        for session_path in sessions_paths:
            session_id = session_path.split(sep)[-1]
            row_df = pd.DataFrame(
                [[participant_id, session_id]], columns=["participant_id", "session_id"]
            )

            # Check t1-volume outputs
            if path.exists(path.join(session_path, "t1", "spm", "segmentation")):
                row_df.loc[0, "t1-volume-segmentation"] = "1"
                for group in groups:
                    group_id = group.split("-")[-1]
                    if path.exists(
                        path.join(
                            session_path,
                            "t1",
                            "spm",
                            "dartel",
                            group,
                            f"{participant_id}_{session_id}_T1w_target-{group_id}_transformation-forward_deformation.nii.gz",
                        )
                    ):
                        row_df.loc[0, f"t1-volume-register-dartel_{group}"] = "1"
                        dartel2mni = glob(
                            path.join(
                                session_path,
                                "t1",
                                "spm",
                                "dartel",
                                group,
                                f"{participant_id}_{session_id}_T1w_segm-*_space-Ixi549Space_modulated-*_probability.nii.gz",
                            )
                        )
                        if len(dartel2mni) > 0:
                            row_df.loc[0, f"t1-volume-dartel2mni_{group}"] = "1"
                            if path.exists(
                                path.join(
                                    session_path,
                                    "t1",
                                    "spm",
                                    "dartel",
                                    group,
                                    "atlas_statistics",
                                )
                            ):
                                row_df.loc[0, f"t1-volume-parcellation_{group}"] = "1"
                            else:
                                row_df.loc[0, f"t1-volume-parcellation_{group}"] = "0"
                        else:
                            row_df.loc[0, f"t1-volume-dartel2mni_{group}"] = "0"
                            row_df.loc[0, f"t1-volume-parcellation_{group}"] = "0"
                    else:
                        row_df.loc[0, f"t1-volume-register-dartel_{group}"] = "0"
                        row_df.loc[0, f"t1-volume-dartel2mni_{group}"] = "0"
                        row_df.loc[0, f"t1-volume-parcellation_{group}"] = "0"
            else:
                row_df.loc[0, "t1-volume-segmentation"] = "0"
                for group in groups:
                    row_df.loc[0, f"t1-volume-register-dartel_{group}"] = "0"
                    row_df.loc[0, f"t1-volume-dartel2mni_{group}"] = "0"
                    row_df.loc[0, f"t1-volume-parcellation_{group}"] = "0"

            # Check t1-linear outputs
            if path.exists(path.join(session_path, "t1_linear")):
                row_df.loc[0, "t1-linear"] = "1"
            else:
                row_df.loc[0, "t1-linear"] = "0"

            # Check t1-freesurfer outputs
            if path.exists(path.join(session_path, "t1", "freesurfer_cross_sectional")):
                row_df.loc[0, "t1-freesurfer"] = "1"
            else:
                row_df.loc[0, "t1-freesurfer"] = "0"

            # Check pet-volume outputs
            for group in groups:
                for trc in trc_avail:
                    for pvc in [True, False]:
                        pet_pattern = path.join(
                            session_path, "pet", "preprocessing", group, f"*{trc}*"
                        )
                        pet_paths = glob(pet_pattern)
                        if pvc:
                            pet_paths = [
                                pet_path for pet_path in pet_paths if "pvc" in pet_path
                            ]
                        else:
                            pet_paths = [
                                pet_path
                                for pet_path in pet_paths
                                if "pvc" not in pet_path
                            ]

                        if len(pet_paths) > 0:
                            row_df.loc[0, f"pet-volume_{trc}_{group}_pvc-{pvc}"] = "1"
                        else:
                            row_df.loc[0, f"pet-volume_{trc}_{group}_pvc-{pvc}"] = "0"

            # Check pet-surface outputs
            for trc in trc_avail:
                pet_pattern = path.join(session_path, "pet", "surface", f"*{trc}*")
                if len(glob(pet_pattern)) > 0:
                    row_df.loc[0, f"pet-surface_{trc}"] = "1"
                else:
                    row_df.loc[0, f"pet-surface_{trc}"] = "0"

            output_df = pd.concat([output_df, row_df])

    output_df.sort_values(["participant_id", "session_id"], inplace=True)
    output_df.to_csv(out_file, sep="\t", index=False)


def compute_missing_mods(bids_dir, out_dir, output_prefix=""):
    """Compute the list of missing modalities for each subject in a BIDS compliant dataset.

    Args:
        bids_dir: path to the BIDS directory
        out_dir: path to the output folder
        output_prefix: string that replace the default prefix ('missing_mods_') in the name of all the output files
    created
    """
    import os
    from glob import glob
    from os import path

    import pandas as pd

    from ..converter_utils import (
        MissingModsTracker,
        print_longitudinal_analysis,
        print_statistics,
    )

    os.makedirs(out_dir, exist_ok=True)

    # Find all the modalities and sessions available for the input dataset
    mods_and_sess = find_mods_and_sess(bids_dir)
    sessions_found = mods_and_sess["sessions"]
    mods_and_sess.pop("sessions")
    mods_avail_dict = mods_and_sess
    mods_avail = [j for i in mods_avail_dict.values() for j in i]
    cols_dataframe = mods_avail[:]
    cols_dataframe.insert(0, "participant_id")
    mmt = MissingModsTracker(sessions_found, mods_avail)

    if output_prefix == "":
        out_file_name = "missing_mods_"
    else:
        out_file_name = output_prefix + "_"

    summary_file = open(path.join(out_dir, out_file_name + "summary.txt"), "w")
    analysis_file = open(path.join(out_dir, "analysis.txt"), "w")
    missing_mods_df = pd.DataFrame(columns=cols_dataframe)
    row_to_append_df = pd.DataFrame(columns=cols_dataframe)
    subjects_paths_lists = glob(path.join(bids_dir, "*sub-*"))
    subjects_paths_lists.sort()

    if len(subjects_paths_lists) == 0:
        raise IOError("No subjects found or dataset not BIDS complaint.")
    # Check the modalities available for each session
    for ses in sessions_found:
        for sub_path in subjects_paths_lists:
            mods_avail_bids = []
            subj_id = sub_path.split(os.sep)[-1]
            row_to_append_df["participant_id"] = pd.Series(subj_id)
            ses_path_avail = glob(path.join(sub_path, ses))
            if len(ses_path_avail) == 0:
                mmt.increase_missing_ses(ses)
                for mod in mods_avail:
                    row_to_append_df[mod] = pd.Series("0")
            else:

                ses_path = ses_path_avail[0]
                mods_paths_folders = glob(path.join(ses_path, "*/"))

                for p in mods_paths_folders:
                    p = p[:-1]
                    mods_avail_bids.append(p.split("/").pop())

                # Check if a modality folder is available and if is empty
                if "func" in mods_avail_bids:
                    # Extract all the task available
                    for m in mods_avail_dict["func"]:
                        tokens = m.split("_")
                        task_name = tokens[1]
                        task_avail_list = glob(
                            path.join(ses_path, "func", "*" + task_name + "*")
                        )

                        if len(task_avail_list) == 0:
                            row_to_append_df[m] = pd.Series("0")
                        else:
                            row_to_append_df[m] = pd.Series("1")
                # If the folder is not available but the modality is
                # in the list of the available one mark it as missing
                else:
                    if "func" in mods_avail_dict:
                        for m in mods_avail_dict["func"]:
                            row_to_append_df[m] = pd.Series("0")
                        mmt.add_missing_mod(ses, m)

                if "dwi" in mods_avail_bids:
                    row_to_append_df["dwi"] = pd.Series("1")
                else:
                    if "dwi" in mods_avail:
                        row_to_append_df["dwi"] = pd.Series("0")
                        mmt.add_missing_mod(ses, "dwi")

                if "anat" in mods_avail_bids:
                    for m in mods_avail_dict["anat"]:
                        anat_aval_list = glob(path.join(ses_path, "anat", "*.nii.gz"))
                        anat_aval_list = [
                            elem for elem in anat_aval_list if m.lower() in elem.lower()
                        ]
                        if len(anat_aval_list) > 0:
                            row_to_append_df[m] = pd.Series("1")
                        else:
                            row_to_append_df[m] = pd.Series("0")
                            mmt.add_missing_mod(ses, m)
                else:
                    if "anat" in mods_avail_dict:
                        for m in mods_avail_dict["anat"]:
                            row_to_append_df[m] = pd.Series("0")
                            mmt.add_missing_mod(ses, m)

                if "fmap" in mods_avail_bids:
                    row_to_append_df["fmap"] = pd.Series("1")
                else:
                    if "fmap" in mods_avail:
                        row_to_append_df["fmap"] = pd.Series("0")
                        mmt.add_missing_mod(ses, "fmap")
                if "pet" in mods_avail_bids:
                    # Extract all the task available
                    for m in mods_avail_dict["pet"]:
                        tokens = m.split("_")
                        pet_acq = tokens[1]
                        acq_avail_list = glob(
                            path.join(ses_path, "pet", "*" + pet_acq + "*")
                        )

                        if len(acq_avail_list) == 0:
                            row_to_append_df[m] = pd.Series("0")
                        else:
                            row_to_append_df[m] = pd.Series("1")
                # If the folder is not available but the modality is
                # in the list of the available one mark it as missing
                else:
                    if "pet" in mods_avail_dict:
                        for m in mods_avail_dict["pet"]:
                            row_to_append_df[m] = pd.Series("0")
                        mmt.add_missing_mod(ses, m)

            missing_mods_df = missing_mods_df.append(row_to_append_df)
            row_to_append_df = pd.DataFrame(columns=cols_dataframe)

        missing_mods_df = missing_mods_df[cols_dataframe]
        missing_mods_df.to_csv(
            path.join(out_dir, out_file_name + ses + ".tsv"),
            sep="\t",
            index=False,
            encoding="utf-8",
        )
        missing_mods_df = pd.DataFrame(columns=cols_dataframe)

    print_statistics(summary_file, len(subjects_paths_lists), sessions_found, mmt)
    print_longitudinal_analysis(
        analysis_file, bids_dir, out_dir, sessions_found, out_file_name
    )


def create_subs_sess_list(
    input_dir, output_dir, file_name=None, is_bids_dir=True, use_session_tsv=False
):
    """Create the file subject_session_list.tsv that contains the list of the visits for each subject for a BIDS or CAPS compliant dataset.

    Args:
        input_dir (str): Path to the BIDS or CAPS directory.
        output_dir (str): Path to the output directory
        file_name: name of the output file
        is_bids_dir (boolean): Specify if input_dir is a BIDS directory or
            not (i.e. a CAPS directory)
        use_session_tsv (boolean): Specify if the list uses the sessions listed in the sessions.tsv files
    """
    import os
    from glob import glob
    from os import path

    import pandas as pd

    os.makedirs(output_dir, exist_ok=True)

    if not file_name:
        file_name = "subjects_sessions_list.tsv"
    subjs_sess_tsv = open(path.join(output_dir, file_name), "w")
    subjs_sess_tsv.write("participant_id" + "\t" + "session_id" + "\n")

    if is_bids_dir:
        path_to_search = input_dir
    else:
        path_to_search = path.join(input_dir, "subjects")
    subjects_paths = glob(path.join(path_to_search, "*sub-*"))

    # Sort the subjects list
    subjects_paths.sort()

    if len(subjects_paths) == 0:
        raise IOError("Dataset empty or not BIDS/CAPS compliant.")

    for sub_path in subjects_paths:
        subj_id = sub_path.split(os.sep)[-1]

        if use_session_tsv:
            session_df = pd.read_csv(
                path.join(sub_path, subj_id + "_sessions.tsv"), sep="\t"
            )
            session_list = list(session_df["session_id"].to_numpy())
            for session in session_list:
                subjs_sess_tsv.write(subj_id + "\t" + session + "\n")

        else:
            sess_list = glob(path.join(sub_path, "*ses-*"))

            for ses_path in sess_list:
                session_name = ses_path.split(os.sep)[-1]
                subjs_sess_tsv.write(subj_id + "\t" + session_name + "\n")

    subjs_sess_tsv.close()


def center_nifti_origin(input_image, output_image):
    """Put the origin of the coordinate system at the center of the image.

    Args:
        input_image: path to the input image
        output_image: path to the output image (where the result will be stored)

    Returns:
        path of the output image created
    """
    import os
    from os.path import isfile

    import nibabel as nib
    import numpy as np
    from nibabel.spatialimages import ImageFileError

    error_str = None
    try:
        img = nib.load(input_image)
    except FileNotFoundError:
        error_str = f"No such file {input_image}"
    except ImageFileError:
        error_str = f"File {input_image} could not be read"

    except Exception as e:
        error_str = f"File {input_image} could not be loaded with nibabel: {e}"

    if not error_str:
        try:
            canonical_img = nib.as_closest_canonical(img)
            hd = canonical_img.header

            qform = np.zeros((4, 4))
            for i in range(1, 4):
                qform[i - 1, i - 1] = hd["pixdim"][i]
                qform[i - 1, 3] = -1.0 * hd["pixdim"][i] * hd["dim"][i] / 2.0
            new_img = nib.Nifti1Image(
                canonical_img.get_data(caching="unchanged"), affine=qform, header=hd
            )

            # Without deleting already-existing file, nib.save causes a severe bug on Linux system
            if isfile(output_image):
                os.remove(output_image)

            nib.save(new_img, output_image)
            if not isfile(output_image):
                error_str = (
                    f"NIfTI file created but Clinica could not save it to {output_image}. "
                    "Please check that the output folder has the correct permissions."
                )
        except Exception as e:
            error_str = (
                "File "
                + input_image
                + " could not be processed with nibabel: "
                + str(e)
            )

    return output_image, error_str


def center_all_nifti(bids_dir, output_dir, modality, center_all_files=False):
    """Center all the NIfTI images of the input BIDS folder into the empty output_dir specified in argument.

    All the files from bids_dir are copied into output_dir, then all the NIfTI images we can found are replaced by their
    centered version if their center if off the origin by more than 50 mm.

    Args:
        bids_dir: (str) path to bids directory
        output_dir: (str) path to EMPTY output directory
        modality: (list of str) modalities to convert
        center_all_files: (bool) center only files that may cause problem for SPM if false. If true, center all NIfTI

    Returns:
        List of the centered files
    """
    from glob import glob
    from os import listdir
    from os.path import basename, isdir, isfile, join
    from shutil import copy, copy2, copytree

    from clinica.utils.exceptions import ClinicaBIDSError
    from clinica.utils.inputs import check_bids_folder

    # output and input must be different, so that we do not mess with user's data
    if bids_dir == output_dir:
        raise ClinicaBIDSError("Input BIDS and output directories must be different")

    # check that input is a BIDS dir
    check_bids_folder(bids_dir)

    for f in listdir(bids_dir):
        if isdir(join(bids_dir, f)) and not isdir(join(output_dir, f)):
            copytree(join(bids_dir, f), join(output_dir, f), copy_function=copy)
        elif isfile(join(bids_dir, f)) and not isfile(join(output_dir, f)):
            copy(join(bids_dir, f), output_dir)

    pattern = join(output_dir, "**/*.nii*")
    nifti_files = glob(pattern, recursive=True)

    # Now filter this list by elements in modality list
    #   For each file:
    #       if any modality name (lowercase) is found in the basename of the file:
    #           keep the file
    nifti_files_filtered = [
        f
        for f in nifti_files
        if any(elem.lower() in basename(f).lower() for elem in modality)
    ]

    # Remove those who are centered
    if not center_all_files:
        nifti_files_filtered = [
            file for file in nifti_files_filtered if not is_centered(file)
        ]

    all_errors = []
    for f in nifti_files_filtered:
        print("Handling " + f)
        _, current_error = center_nifti_origin(f, f)
        if current_error:
            all_errors.append(current_error)
    if len(all_errors) > 0:
        final_error_msg = (
            "Clinica encoutered "
            + str(len(all_errors))
            + " error(s) while trying to center all NIfTI images.\n"
        )
        for error in all_errors:
            final_error_msg += "\n" + error
        raise RuntimeError(final_error_msg)
    return nifti_files_filtered


def are_far_appart(file1, file2, threshold=80):
    """Tell if 2 files have a center located at more than a threshold distance.

    Args:
        file1: (str) path to the first nifti file
        file2: (str) path to the second nifti file
        threshold: threshold to consider whether 2 files are too far appart

    Returns:
        True if distance between `file1` and `file2` is greter than `threshold`, False otherwise.
    """
    from os.path import isfile

    import numpy as np

    assert isfile(file1)
    assert isfile(file2)

    center1 = get_world_coordinate_of_center(file1)
    center2 = get_world_coordinate_of_center(file2)

    return np.linalg.norm(center2 - center1, ord=2) > threshold


def write_list_of_files(file_list, output_file):
    """Save `file_list` list of files into `output_file` text file.

    Args:
        file_list: (list of str) of path to files
        output_file: (str) path to the output txt file

    Returns:
        output_file
    """
    from os.path import isfile

    assert isinstance(file_list, list), "First argument must be a list"
    assert isinstance(output_file, str), "Second argument must be a str"
    if isfile(output_file):
        return None

    text_file = open(output_file, "w+")
    for created_file in file_list:
        text_file.write(created_file + "\n")
    text_file.close()
    return output_file


def check_relative_volume_location_in_world_coordinate_system(
    label_1, nifti_list1, label_2, nifti_list2, bids_dir, modality
):
    """
    Check if the NIfTI file list nifti_list1 and nifti_list2 provided in argument are not too far apart (otherwise coreg
    in SPM may fail. Norm between center of volumes of 2 files must be less than 80 mm.

    Args:
        label_1: label of the first nifti_list1 files (used in potential warning message)
        nifti_list1: first set of files
        label_2: label of the second nifti_list
        nifti_list2: second set of files, must be same length as nifti_list1
        bids_dir: bids directory (used in potential warning message)
        modality: string that must be used in argument of: clinica iotools bids --modality MODALITY (used in potential
                warning message)
    """
    import sys
    from os.path import abspath, basename

    import numpy as np

    from clinica.utils.stream import cprint

    center_coordinate_1 = [get_world_coordinate_of_center(file) for file in nifti_list1]
    center_coordinate_2 = [get_world_coordinate_of_center(file) for file in nifti_list2]

    l2_norm = [
        np.linalg.norm(center_1 - center_2)
        for center_1, center_2 in zip(center_coordinate_1, center_coordinate_2)
    ]
    pairs_with_problems = [i for i, norm in enumerate(l2_norm) if norm > 80]

    if len(pairs_with_problems) > 0:
        warning_message = (
            f"It appears that {str(len(pairs_with_problems))} "
            "pairs of files have an important relative offset. "
            "SPM coregistration has a high probability to fail on these files:\n\n"
        )

        # File column width : 3 spaces more than the longest string to display
        file1_width = max(
            3 + len(label_1),
            3
            + max(
                len(basename(file))
                for file in [nifti_list1[k] for k in pairs_with_problems]
            ),
        )
        file2_width = max(
            3 + len(label_2),
            3
            + max(
                len(basename(file))
                for file in [nifti_list2[k] for k in pairs_with_problems]
            ),
        )

        norm_width = len("Relative distance")

        warning_message += (
            "%-"
            + str(file1_width)
            + "s%-"
            + str(file2_width)
            + "s%-"
            + str(norm_width)
            + "s"
        ) % (label_1, label_2, "Relative distance")

        warning_message += "\n" + "-" * (file1_width + file2_width + norm_width) + "\n"
        for file1, file2, norm in zip(
            [nifti_list1[k] for k in pairs_with_problems],
            [nifti_list2[k] for k in pairs_with_problems],
            [l2_norm[k] for k in pairs_with_problems],
        ):
            # Nice formatting as array
            # % escape character
            # - aligned to the left, with the size of the column
            # s = string, f = float
            # . for precision with float
            # https://docs.python.org/2/library/stdtypes.html#string-formatting for more information
            warning_message += (
                "%-"
                + str(file1_width)
                + "s%-"
                + str(file2_width)
                + "s%-"
                + str(norm_width)
                + ".2f\n"
            ) % (str(basename(file1)), str(basename(file2)), norm)
        warning_message += (
            "\nClinica provides a tool to counter this problem by replacing the center "
            "of the volume at the origin of the world coordinates.\nUse the following "
            "command line to correct the header of the faulty NIFTI volumes in a new folder:\n\n"
            f"`clinica iotools center-nifti {abspath(bids_dir)} {abspath(bids_dir)}_centered --modality {modality}`\n\n"
            "You will find more information on the command by typing `clinica iotools center-nifti` in the console."
        )
        cprint(msg=warning_message, lvl="warning")
        if not click.confirm("Do you still want to launch the pipeline?"):
            click.echo("Clinica will now exit...")
            sys.exit(0)


def check_volume_location_in_world_coordinate_system(
    nifti_list, bids_dir, modality="t1w", skip_question=False
):
    """
    Check if the NIfTI file list nifti_list provided in argument are aproximately centered around the origin of the
    world coordinates. (Problem may arise with SPM segmentation

    If yes, we warn the user of this problem, and propose him to exit clinica in order for him to run:
        clinica iotools center-nifti ...
    or to continue with the execution of the pipeline

    Args:
        nifti_list: (list of str) list of path to nifti files
        bids_dir: (str) path to bids directory associated with this check (in order to propose directly the good
            command line for center-nifti tool)
        modality: (str) to propose directly the good command line option
        skip_question: (bool) if True user input is not asked for and the answer is automatically yes
    """
    import sys
    from os.path import abspath, basename

    import click
    import numpy as np

    list_non_centered_files = [file for file in nifti_list if not is_centered(file)]
    if len(list_non_centered_files) > 0:
        centers = [
            get_world_coordinate_of_center(file) for file in list_non_centered_files
        ]
        l2_norm = [np.linalg.norm(center, ord=2) for center in centers]

        # File column width : 3 spaces more than the longest string to display
        file_width = 3 + max(len(basename(file)) for file in list_non_centered_files)
        # Center column width (with a fixed minimum size) : 3 spaces more than the longest string to display
        center_width = max(
            len("Coordinate of center") + 3,
            3 + max(len(str(center)) for center in centers),
        )

        warning_message = (
            f"It appears that {str(len(list_non_centered_files))} files "
            "have a center way out of the origin of the world coordinate system. SPM has a high "
            "probability to fail on these files (for coregistration or segmentation):\n\n"
        )
        warning_message += (
            "%-" + str(file_width) + "s%-" + str(center_width) + "s%-s"
        ) % ("File", "Coordinate of center", "Distance to origin")
        # 18 is the length of the string 'Distance to origin'
        warning_message += "\n" + "-" * (file_width + center_width + 18) + "\n"
        for file, center, l2 in zip(list_non_centered_files, centers, l2_norm):
            # Nice formatting as array
            # % escape character
            # - aligned to the left, with the size of the column
            # s = string, f = float
            # . for precision with float
            # https://docs.python.org/2/library/stdtypes.html#string-formatting for more information
            warning_message += (
                "%-" + str(file_width) + "s%-" + str(center_width) + "s%-25.2f\n"
            ) % (basename(file), str(center), l2)

        cmd_line = f"`clinica iotools center-nifti {abspath(bids_dir)} {abspath(bids_dir)}_centered --modality {modality}`"

        warning_message += (
            "\nIf you are trying to launch the t1-freesurfer pipeline, you can ignore this message "
            "if you do not want to run the pet-surface pipeline afterward."
        )

        warning_message += (
            "\nClinica provides a tool to counter this problem by replacing the center of the volume"
            " at the origin of the world coordinates.\nUse the following command line to correct the "
            f"header of the faulty NIFTI volumes in a new folder:\n{cmd_line}"
            "You will find more information on the command by typing "
            "clinica iotools center-nifti in the console."
        )

        click.echo(warning_message)

        if not skip_question:
            if not click.confirm("Do you still want to launch the pipeline?"):
                click.echo("Clinica will now exit...")
                sys.exit(0)


def is_centered(nii_volume, threshold_l2=50):
    """Tell if a NIfTI volume is centered on the origin of the world coordinate system.

    SPM has troubles to segment files if the center of the volume is not close from the origin of the world coordinate
    system. A series of experiment have been conducted: we take a volume whose center is on the origin of the world
    coordinate system. We add an offset using coordinates of affine matrix [0, 3], [1, 3], [2, 3] (or by modifying the
    header['srow_x'][3], header['srow_y'][3], header['srow_z'][3], this is strictly equivalent).

    It has been determined that volumes were still segmented with SPM when the L2 distance between origin and center of
    the volume did not exceed 100 mm. Above this distance, either the volume is either not segmented (SPM error), or the
    produced segmentation is wrong (not the shape of a brain anymore)

    Args:
        nii_volume: path to NIfTI volume
        threshold_l2: maximum distance between origin of the world coordinate system and the center of the volume to
                    be considered centered. The threshold were SPM segmentation stops working is around 100 mm
                    (it was determined empirically after several trials on a genrated dataset), so default value is 50
                    mm in order to have a security margin, even when dealing with coregistred files afterward)

    Returns:
        True or False
    """
    import numpy as np

    center = get_world_coordinate_of_center(nii_volume)

    # Compare to the threshold and retun boolean
    # if center is a np.nan, comparison will be False, and False will be returned
    distance_from_origin = np.linalg.norm(center, ord=2)
    # if not np.isnan(distance_from_origin):
    #     print('\t' + basename(nii_volume) + ' has its center at {0:.2f} mm of the origin.'.format(distance_from_origin))
    if distance_from_origin < threshold_l2:
        return True
    else:
        # If center is a np.nan,
        return False


def get_world_coordinate_of_center(nii_volume):
    """Extract the world coordinates of the center of the image.

    Based on methods described here: https://brainder.org/2012/09/23/the-nifti-file-format/

    Args:
        nii_volume: path to nii volume

    Returns:
        [Returns]
    """
    from os.path import isfile

    import nibabel as nib
    import numpy as np

    from clinica.utils.stream import cprint

    assert isinstance(nii_volume, str), "input argument nii_volume must be a str"
    assert isfile(nii_volume), "input argument must be a path to a file"

    try:
        orig_nifti = nib.load(nii_volume)
    except nib.filebasedimages.ImageFileError:
        cprint(
            msg=f"File {nii_volume} could not be read by nibabel. Is it a valid NIfTI file ?",
            lvl="warning",
        )
        return np.nan

    head = orig_nifti.header

    if isinstance(head, nib.freesurfer.mghformat.MGHHeader):
        # If MGH volume
        center_coordinates_world = vox_to_world_space_method_3_bis(
            head["dims"][0:3] / 2, head
        )
    else:
        # Standard NIfTI volume
        center_coordinates = get_center_volume(head)

        if head["qform_code"] > 0:
            center_coordinates_world = vox_to_world_space_method_2(
                center_coordinates, head
            )
        elif head["sform_code"] > 0:
            center_coordinates_world = vox_to_world_space_method_3(
                center_coordinates, head
            )
        elif head["sform_code"] == 0:
            center_coordinates_world = vox_to_world_space_method_1(
                center_coordinates, head
            )
        else:
            center_coordinates_world = np.nan
    return center_coordinates_world


def get_center_volume(header):
    """Get the voxel coordinates of the center of the data, using header information.

    Args:
        header: a nifti header

    Returns:
        Voxel coordinates of the center of the volume
    """
    import numpy as np

    center_x = header["dim"][1] / 2
    center_y = header["dim"][2] / 2
    center_z = header["dim"][3] / 2
    return np.array([center_x, center_y, center_z])


def vox_to_world_space_method_1(coordinates_vol, header):
    """
    The Method 1 is for compatibility with analyze and is not supposed to be used as the main orientation method. But it
    is used if sform_code = 0. The world coordinates are determined simply by scaling by the voxel size by their
    dimension stored in pixdim. More information here: https://brainder.org/2012/09/23/the-nifti-file-format/
    Args:
        coordinates_vol: coordinate in the volume (raw data)
        header: header object

    Returns:
        Coordinates in the world space
    """
    import numpy as np

    return np.array(coordinates_vol) * np.array(
        header["pixdim"][1], header["pixdim"][2], header["pixdim"][3]
    )


def vox_to_world_space_method_2(coordinates_vol, header):
    """
    The Method 2 is used when short qform_code is larger than zero. To get the coordinates, we multiply a rotation
    matrix (r_mat) by coordinates_vol, then perform hadamart with pixel dimension pixdim (like in method 1). Then we add
    an offset (qoffset_x, qoffset_y, qoffset_z)

    Args:
        coordinates_vol: coordinate in the volume (raw data)
        header: header object

    Returns:
        Coordinates in the world space
    """
    import numpy as np

    def get_r_matrix(h):
        """Get rotation matrix.

        More information here: https://brainder.org/2012/09/23/the-nifti-file-format/

        Args:
            h: header

        Returns:
            Rotation matrix
        """
        b = h["quatern_b"]
        c = h["quatern_c"]
        d = h["quatern_d"]
        a = np.sqrt(1 - (b ** 2) - (c ** 2) - (d ** 2))
        r = np.zeros((3, 3))
        r[0, 0] = (a ** 2) + (b ** 2) - (c ** 2) - (d ** 2)
        r[0, 1] = 2 * ((b * c) - (a * d))
        r[0, 2] = 2 * ((b * d) + (a * c))
        r[1, 0] = 2 * ((b * c) + (a * d))
        r[1, 1] = (a ** 2) + (c ** 2) - (b ** 2) - (d ** 2)
        r[1, 2] = 2 * ((c * d) - (a * b))
        r[2, 0] = 2 * ((b * d) - (a * c))
        r[2, 1] = 2 * ((b * d) - (a * c))
        r[2, 2] = (a ** 2) + (d ** 2) - (b ** 2) - (c ** 2)
        return r

    i = coordinates_vol[0]
    j = coordinates_vol[1]
    k = coordinates_vol[2]
    if header["qform_code"] > 0:
        r_mat = get_r_matrix(header)
    else:
        # Should never be reached
        raise ValueError("qform_code must be greater than 0 to use this method")
    q = header["pixdim"][0]
    if q not in [-1, 1]:
        print("q was " + str(q), ", now is 1")
        q = 1
    return np.dot(r_mat, np.array([i, j, q * k])) * np.array(
        header["pixdim"][1:4]
    ) + np.array([header["qoffset_x"], header["qoffset_y"], header["qoffset_z"]])


def vox_to_world_space_method_3(coordinates_vol, header):
    """
    This method is used when sform_code is larger than zero. It relies on a full affine matrix, stored in the header in
    the fields srow_[x,y,y], to map voxel to world coordinates.
    When a nifti file is created with raw data and affine=..., this is this method that is used to decypher the
    voxel-to-world correspondance.

    Args:
        coordinates_vol: coordinate in the volume (raw data)
        header: header object

    Returns:
        Coordinates in the world space
    """
    import numpy as np

    def get_aff_matrix(h):
        """Get affine transformation matrix.

        See details here: https://brainder.org/2012/09/23/the-nifti-file-format/

        Args:
            h: header

        Returns:
            affine transformation matrix
        """
        mat = np.zeros((4, 4))
        mat[0, 0] = h["srow_x"][0]
        mat[0, 1] = h["srow_x"][1]
        mat[0, 2] = h["srow_x"][2]
        mat[0, 3] = h["srow_x"][3]
        mat[1, 0] = h["srow_y"][0]
        mat[1, 1] = h["srow_y"][1]
        mat[1, 2] = h["srow_y"][2]
        mat[1, 3] = h["srow_y"][3]
        mat[2, 0] = h["srow_z"][0]
        mat[2, 1] = h["srow_z"][1]
        mat[2, 2] = h["srow_z"][2]
        mat[2, 3] = h["srow_z"][3]
        mat[3, 3] = 1
        return mat

    if header["sform_code"] > 0:
        aff = get_aff_matrix(header)
    else:
        # Should never be reached
        raise ValueError("sform_code has a value > 0, so method 3 cannot be used")

    homogeneous_coord = np.concatenate(
        (np.array(coordinates_vol), np.array([1])), axis=0
    )
    return np.dot(aff, homogeneous_coord)[0:3]


def vox_to_world_space_method_3_bis(coordinates_vol, header):
    """
    This method relies on the same technique as method 3, but for images created by FreeSurfer (MGHImage, MGHHeader).
    Args:
        coordinates_vol: coordinate in the volume (raw data)
        header: nib.freesurfer.mghformat.MGHHeader object

    Returns:
        Coordinates in the world space
    """
    import numpy as np

    affine_trensformation_matrix = header.get_affine()
    homogeneous_coord = np.concatenate(
        (np.array(coordinates_vol), np.array([1])), axis=0
    )
    return np.dot(affine_trensformation_matrix, homogeneous_coord)[0:3]
