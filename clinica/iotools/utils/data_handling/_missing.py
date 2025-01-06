from dataclasses import dataclass
from os import PathLike
from pathlib import Path
from typing import List

import pandas as pd

__all__ = [
    "compute_missing_mods",
    "compute_missing_processing",
]


@dataclass
class SubjectSession:
    participant_id: str
    session_id: str
    session_path: Path

    @property
    def sub_ses_id(self) -> str:
        return f"{self.participant_id}_{self.session_id}"


def compute_missing_mods(
    bids_dir: PathLike, out_dir: PathLike, output_prefix: str = ""
) -> None:
    """Compute the list of missing modalities for each subject in a BIDS compliant dataset.

    Parameters
    ----------
    bids_dir : PathLike
        Path to the BIDS directory.

    out_dir : PathLike
        Path to the output folder.

    output_prefix : str, optional
        String that replaces the default prefix ('missing_mods_')
        in the name of all the created output files.
        Default = "".
    """
    import os
    from glob import glob
    from os import path
    from pathlib import Path

    import pandas as pd

    from clinica.iotools.converter_utils import (
        MissingModsTracker,
        write_longitudinal_analysis,
        write_statistics,
    )

    out_dir = Path(out_dir)
    bids_dir = Path(bids_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Find all the modalities and sessions available for the input dataset
    mods_and_sess = _find_mods_and_sess(bids_dir)
    sessions_found = mods_and_sess["sessions"]
    mods_and_sess.pop("sessions")
    mods_avail_dict = mods_and_sess
    mods_avail = [j for i in mods_avail_dict.values() for j in i]
    cols_dataframe = mods_avail[:]
    cols_dataframe.insert(0, "participant_id")
    mmt = MissingModsTracker(sessions_found, mods_avail)

    out_file_name = "missing_mods_" if output_prefix == "" else output_prefix + "_"

    missing_mods_df = pd.DataFrame(columns=cols_dataframe)
    row_to_append_df = pd.DataFrame(columns=cols_dataframe)
    subjects_paths_lists = glob(path.join(bids_dir, "*sub-*"))
    subjects_paths_lists.sort()

    if len(subjects_paths_lists) == 0:
        raise IOError("No subjects found or dataset not BIDS compliant.")
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

            missing_mods_df = pd.concat([missing_mods_df, row_to_append_df])
            row_to_append_df = pd.DataFrame(columns=cols_dataframe)

        missing_mods_df = missing_mods_df[cols_dataframe]
        missing_mods_df.to_csv(
            path.join(out_dir, out_file_name + ses + ".tsv"),
            sep="\t",
            index=False,
            encoding="utf-8",
        )
        missing_mods_df = pd.DataFrame(columns=cols_dataframe)

    write_statistics(
        out_dir / (out_file_name + "summary.txt"),
        len(subjects_paths_lists),
        sessions_found,
        mmt,
    )
    write_longitudinal_analysis(
        out_dir / "analysis.txt", bids_dir, out_dir, sessions_found, out_file_name
    )


def compute_missing_processing(
    bids_dir: PathLike, caps_dir: PathLike, out_file: PathLike
):
    """Compute the list of missing processing for each subject in a CAPS compliant dataset.

    Parameters
    ----------
    bids_dir : PathLike
        Path to the BIDS directory.

    caps_dir : PathLike
        Path to the CAPS directory.

    out_file : PathLike
        Path to the output file (filename included).
    """
    bids_dir = Path(bids_dir)
    caps_dir = Path(caps_dir)
    output_df = pd.DataFrame()
    groups = _get_groups(caps_dir)
    tracers = _get_pet_tracers(bids_dir)

    for subject_path in (caps_dir / "subjects").glob("sub-*"):
        participant_id = subject_path.parent.name
        for session_path in subject_path.glob("ses-*"):
            session_id = session_path.parent.name
            subject = SubjectSession(participant_id, session_id, session_path)
            row = _compute_missing_processing_single_row(subject, groups, tracers)
            output_df = pd.concat([output_df, row])
    output_df.sort_values(["participant_id", "session_id"], inplace=True)
    output_df.to_csv(out_file, sep="\t", index=False)


def _get_groups(caps_dir: Path) -> List[str]:
    if (caps_dir / "groups").exists():
        return [f.name for f in (caps_dir / "groups").iterdir()]
    return []


def _get_pet_tracers(bids_dir: Path) -> List[str]:
    """Retrieve pet tracers available."""
    from clinica.utils.stream import cprint

    modalities = _find_mods_and_sess(bids_dir)
    modalities.pop("sessions")
    tracers = [
        j.split("_")[1]
        for modality in modalities.values()
        for j in modality
        if "pet" in j
    ]
    cprint(f"Available tracers: {tracers}", lvl="info")

    return tracers


def _compute_missing_processing_single_row(
    subject: SubjectSession,
    groups: List[str],
    tracers: List[str],
) -> pd.DataFrame:
    """Compute a single row of the missing processing dataframe."""
    row = pd.DataFrame(
        [[subject.participant_id, subject.session_id]],
        columns=["participant_id", "session_id"],
    )
    row = _compute_missing_processing_t1_volume(row, groups, subject)
    for pipeline in ("t1-linear", "flair_linear"):
        row[0, pipeline] = "1" if (subject.session_path / pipeline).exists() else "0"
    row[0, "t1-freesurfer"] = (
        "1"
        if (subject.session_path / "t1" / "freesurfer_cross_sectional").exists()
        else "0"
    )
    row = _compute_missing_processing_pet_volume(
        row, groups, tracers, subject.session_path
    )
    for tracer in tracers:
        row.loc[0, f"pet-surface_{tracer}"] = (
            "1"
            if len(
                [
                    f
                    for f in (subject.session_path / "pet" / "surface").glob(
                        f"*{tracer}*"
                    )
                ]
            )
            > 0
            else "0"
        )
    for tracer in tracers:
        row.loc[0, f"pet-linear_{tracer}"] = (
            "1"
            if len(
                [f for f in (subject.session_path / "pet_linear").glob(f"*{tracer}*")]
            )
            > 0
            else "0"
        )
    return row


def _compute_missing_processing_t1_volume(
    df: pd.DataFrame,
    groups: List[str],
    subject: SubjectSession,
) -> pd.DataFrame:
    if (subject.session_path / "t1" / "spm" / "segmentation").exists():
        df.loc[0, "t1-volume-segmentation"] = "1"
        for group in groups:
            group_id = group.split("-")[-1]
            filename = f"{subject.sub_ses_id}_T1w_target-{group_id}_transformation-forward_deformation.nii.gz"
            if (
                subject.session_path / "t1" / "spm" / "dartel" / group / filename
            ).exists():
                df.loc[0, f"t1-volume-register-dartel_{group}"] = "1"
                pattern = f"{subject.sub_ses_id}_T1w_segm-*_space-Ixi549Space_modulated-*_probability.nii.gz"
                dartel2mni = [
                    f
                    for f in (
                        subject.session_path / "t1" / "spm" / "dartel" / group
                    ).glob(pattern)
                ]
                if len(dartel2mni) > 0:
                    df.loc[0, f"t1-volume-dartel2mni_{group}"] = "1"
                    if (
                        subject.session_path
                        / "t1"
                        / "spm"
                        / "dartel"
                        / group
                        / "atlas_statistics"
                    ).exists():
                        df.loc[0, f"t1-volume-parcellation_{group}"] = "1"
                    else:
                        df.loc[0, f"t1-volume-parcellation_{group}"] = "0"
                else:
                    df.loc[0, f"t1-volume-dartel2mni_{group}"] = "0"
                    df.loc[0, f"t1-volume-parcellation_{group}"] = "0"
            else:
                df.loc[0, f"t1-volume-register-dartel_{group}"] = "0"
                df.loc[0, f"t1-volume-dartel2mni_{group}"] = "0"
                df.loc[0, f"t1-volume-parcellation_{group}"] = "0"
    else:
        df.loc[0, "t1-volume-segmentation"] = "0"
        for group in groups:
            df.loc[0, f"t1-volume-register-dartel_{group}"] = "0"
            df.loc[0, f"t1-volume-dartel2mni_{group}"] = "0"
            df.loc[0, f"t1-volume-parcellation_{group}"] = "0"
    return df


def _compute_missing_processing_pet_volume(
    df: pd.DataFrame,
    groups: List[str],
    tracers: List[str],
    session_path: Path,
) -> pd.DataFrame:
    for group in groups:
        for tracer in tracers:
            pet_paths_pvc = [
                f
                for f in (session_path / "pet" / "preprocessing" / group).glob(
                    f"*{tracer}*"
                )
                if "pvc" in f
            ]
            pet_paths_no_pvc = [
                f
                for f in (session_path / "pet" / "preprocessing" / group).glob(
                    f"*{tracer}*"
                )
                if "pvc" not in f
            ]
            df.loc[0, f"pet-volume_{tracer}_{group}_pvc-True"] = (
                "1" if len(pet_paths_pvc) > 0 else "0"
            )
            df.loc[0, f"pet-volume_{tracer}_{group}_pvc-False"] = (
                "1" if len(pet_paths_no_pvc) > 0 else "0"
            )
    return df


def _find_mods_and_sess(bids_dir: Path) -> dict:
    """Find all the modalities and sessions available for a given BIDS dataset.

    Parameters
    ----------
    bids_dir : Path
        The path to the BIDS dataset.

    Returns
    -------
    mods_dict : dict
        A dictionary that stores the sessions and modalities found and has the following structure.
        {
            'sessions': ['ses-M000', 'ses-M018'],
            'fmap': ['fmap'],
            'anat': ['flair', 't1w'],
            'func': ['func_task-rest'],
            'dwi': ['dwi']
        }
    """
    from collections import defaultdict

    mods_dict = defaultdict(set)

    for sub_path in bids_dir.glob("*sub-*"):
        for session in sub_path.glob("*ses-*"):
            mods_dict["sessions"].add(session.name)
            for modality in session.glob("*/"):
                if modality.name == "func":
                    for func_path in modality.glob("*bold.nii.gz"):
                        func_task = func_path.name.split("_")[2]
                        mods_dict["func"].add(f"func_{func_task}")
                if modality.name == "dwi":
                    mods_dict["dwi"].add("dwi")
                if modality.name == "fmap":
                    mods_dict["fmap"].add("fmap")
                if modality.name == "pet":
                    for pet_path in modality.glob("*pet.nii.gz"):
                        pet_name = pet_path.name.split(".")[0]
                        pet_acq = pet_name.split("_")[2]
                        mods_dict["pet"].add(f"pet_{pet_acq}")
                if modality.name == "anat":
                    for anat_file in modality.glob("*"):
                        if ".nii.gz" in anat_file.name:
                            anat_name = anat_file.name.replace(".nii.gz", "")
                            anat_ext = "nii.gz"
                        else:
                            anat_name = anat_file.stem
                            anat_ext = anat_file.suffix.lstrip(".")
                        if anat_ext != "json":
                            file_parts = anat_name.split("_")
                            anat_type = file_parts[-1].lower()
                            mods_dict["anat"].add(anat_type)

    return mods_dict
