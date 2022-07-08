# coding: utf8

import os
from os import PathLike
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Tuple, Dict, List


def likeliness_measure(
        file1: PathLike,
        file2: PathLike,
        threshold1: Tuple,
        threshold2: Tuple,
        display: bool = False,
 ) -> bool:
    """Compares 2 Nifti inputs, with 2 different thresholds.

    Parameters
    ----------
    file1: Path to first nifti input
    file2: Path to second nifti to compare
    threshold1: Defines the first criteria to meet: threshold1[0] defines the relative
        difference between 2 voxels to be considered different (ex: 1e-4). threshold[1] defines
        the maximum proportion of voxels that can different for the test to be negative.
    threshold2: Defines the second criteria to meet.
    display: If set to True, will display a useful graph to determine optimal threshold for the
        comparison.

    Returns
    -------
    bool
        True if file1 and file2 can be considered similar enough (meeting criterion expressed in threshold1
        and threshold2). False otherwise.
    """
    import matplotlib.pyplot as plt
    import nibabel as nib

    print(" ** comparing " + os.path.basename(file1) + " **")
    data1 = nib.load(file1).get_fdata(dtype="float32")
    data1[np.isnan(data1)] = 0

    data2 = nib.load(file2).get_fdata(dtype="float32")
    data2[np.isnan(data2)] = 0

    # Get mask where data are 0 in data1 and data2
    mask = (data1 == 0) & (data2 == 0)
    data1[mask] = 1
    data2[mask] = 1
    metric = (2 * np.abs(data1 - data2)) / (np.abs(data1) + np.abs(data2))
    metric_flattened = np.ndarray.flatten(metric)

    # Display fig
    if display:
        thresholds = np.logspace(-8, 0, 20)
        percents = np.array(
            [np.sum((metric_flattened > _)) / metric_flattened.size for _ in thresholds]
        )
        fig, ax = plt.subplots()
        ax.semilogx(thresholds, percents)
        ax.grid()
        plt.xlabel("Threshold of relative difference")
        plt.ylabel("Proportion of different voxels")
        plt.show()

    mask_different_voxels_cond1 = metric_flattened > threshold1[0]
    mask_different_voxels_cond2 = metric_flattened > threshold2[0]
    return (
        np.sum(mask_different_voxels_cond1) / metric_flattened.size < threshold1[1]
    ) & (np.sum(mask_different_voxels_cond2) / metric_flattened.size < threshold2[1])


def similarity_measure(
        file1: PathLike,
        file2: PathLike,
        threshold: float,
) -> bool:
    """Compares 2 Nifti inputs using a correlation metric.

    Nifti are equals if the correlation is higher than the specified threshold.

    Parameters
    ----------
    file1: Path to first nifti input
    file2: Path to second nifti to compare
    threshold: Threshold value to be used in the comparison

    Returns
    -------
    bool
        True if file1 and file2 can be considered similar enough, i.e. the
        correlation is higher than the threshold.
    """
    import nipype.pipeline.engine as npe
    from nipype.algorithms.metrics import Similarity

    # Node similarity (nipy required)
    img_similarity = npe.Node(name="img_similarity", interface=Similarity())
    img_similarity.inputs.metric = "cc"  # stands for correlation coefficient
    img_similarity.inputs.volume1 = file1
    img_similarity.inputs.volume2 = file2
    res = img_similarity.run()

    return np.mean(res.outputs.similarity) > threshold


def identical_subject_list(
        sub_ses_list1: PathLike,
        sub_ses_list2: PathLike,
) -> bool:
    """Ensures that both subject_session files are describing the same list.

    Parameters
    ----------
    sub_ses_list1: Path to first nifti input
    sub_ses_list2: Path to second nifti to compare

    Returns
    -------
    bool
        True if sub_ses_list1 and sub_ses_list2 contains the same sessions.
    """

    def is_included(list1: PathLike, list2: PathLike) -> bool:
        readlist1 = pd.read_csv(list1, sep="\t")
        readlist2 = pd.read_csv(list2, sep="\t")

        # If columns are different, files are different
        if list(readlist1.columns) != list(readlist2.columns):
            return False
        else:

            # Extract subject and corresponding session names
            subjects1 = list(readlist1.participant_id)
            sessions1 = list(readlist1.session_id)
            subjects2 = list(readlist2.participant_id)
            sessions2 = list(readlist2.session_id)

            if len(subjects1) != len(subjects2):
                return False
            else:
                for i in range(len(subjects1)):
                    current_sub = subjects1[i]
                    current_ses = sessions1[i]
                    # Compute all the indices in the second list corresponding to the current subject
                    idx_same_sub = [
                        j for j in range(len(subjects2)) if subjects2[j] == current_sub
                    ]
                    if len(idx_same_sub) == 0:  # Current subject not found in
                        return False
                    ses_same_sub = [sessions2[idx] for idx in idx_same_sub]
                    if current_ses not in ses_same_sub:
                        return False
        return True

    # The operation is performed both sides because is_included(list1, list2) != is_included(list2, list1)
    return is_included(sub_ses_list1, sub_ses_list2) & is_included(
        sub_ses_list2, sub_ses_list1
    )


def same_missing_modality_tsv(file1: PathLike, file2: PathLike) -> bool:
    """Compare 2 TSV files generated by the iotool ComputeMissingModalities.

    Only fields participant_id, pet, t1w, func_task - rest are compared. Line order does not matter.

    Parameters
    ----------
    file1: Path to first tsv
    file2: Path to second tsv

    Returns
    -------
    bool
        True if file1 and file2 contains the same information
    """
    df1 = pd.read_csv(file1, sep="\t")
    df2 = pd.read_csv(file2, sep="\t")

    # Extract data and form lists for both files
    subjects1 = list(df1.participant_id)
    pet_AV45_1 = list(df1["pet_trc-18FAV45"])
    pet_FDG_1 = list(df1["pet_trc-18FFDG"])
    t1w1 = list(df1.t1w)
    func_task_rest1 = list(df1["func_task-rest"])

    subjects2 = list(df2.participant_id)
    pet_AV45_2 = list(df2["pet_trc-18FAV45"])
    pet_FDG_2 = list(df2["pet_trc-18FFDG"])
    t1w2 = list(df2.t1w)
    func_task_rest2 = list(df2["func_task-rest"])

    # Subjects are sorted in alphabetical order. The same permutation of element is applied on each column
    subjects1_sorted, pet_AV45_1 = (
        list(t) for t in zip(*sorted(zip(subjects1, pet_AV45_1)))
    )
    subjects2_sorted, pet_AV45_2 = (
        list(t) for t in zip(*sorted(zip(subjects2, pet_AV45_2)))
    )
    subjects1_sorted, pet_FDG_1 = (
        list(t) for t in zip(*sorted(zip(subjects1, pet_FDG_1)))
    )
    subjects2_sorted, pet_FDG_2 = (
        list(t) for t in zip(*sorted(zip(subjects2, pet_FDG_2)))
    )
    subjects1_sorted, t1w1 = (list(t) for t in zip(*sorted(zip(subjects1, t1w1))))
    subjects2_sorted, t1w2 = (list(t) for t in zip(*sorted(zip(subjects2, t1w2))))
    subjects1_sorted, func_task_rest1 = (
        list(t) for t in zip(*sorted(zip(subjects1, func_task_rest1)))
    )
    subjects2_sorted, func_task_rest2 = (
        list(t) for t in zip(*sorted(zip(subjects2, func_task_rest2)))
    )

    # Test is positive when all the sorted list s are equals
    return (
        (subjects1_sorted == subjects2_sorted)
        & (pet_AV45_1 == pet_AV45_2)
        & (pet_FDG_1 == pet_FDG_2)
        & (t1w1 == t1w2)
        & (func_task_rest1 == func_task_rest2)
    )


def compare_folders(outdir: Path, refdir: Path, tmp_path: Path) -> bool:
    from filecmp import cmp

    file_out = tmp_path / "file_out.txt"
    file_ref = tmp_path / "file_ref.txt"
    tree(outdir, file_out)
    tree(refdir, file_ref)

    if not cmp(file_out, file_ref):
        with open(file_out, "r") as fin:
            out_message = fin.read()
        with open(file_ref, "r") as fin:
            ref_message = fin.read()
        raise ValueError(
            "Comparison of out and ref directories shows mismatch :\n "
            "OUT :\n" + out_message + "\n REF :\n" + ref_message
        )
    return True


def tree(dir: Path, file_out: Path):
    """Creates a file (file_out) with a visual tree representing the file
    hierarchy at a given directory

    .. note::
        Does not display empty directories.

    """
    print(type(dir))
    file_content = ""
    for path in sorted(dir.rglob("*")):
        if path.is_dir() and not any(path.iterdir()):
            continue
        depth = len(path.relative_to(dir).parts)
        spacer = "    " * depth
        file_content = file_content + f"{spacer}+ {path.name}\n"
    print(file_content)
    file_out.write_text(file_content)


def clean_folder(path: PathLike, recreate: bool = True):
    """Clean folders under provided path.

    Parameters
    ----------
    path: Path of folder to clean.
    recreate: Whether or not to restore folder
        structure after cleaning.
    """
    from os.path import abspath, exists
    from shutil import rmtree

    abs_path = abspath(path)
    if exists(abs_path):
        rmtree(abs_path)
    if recreate:
        os.makedirs(abs_path)


def list_files_with_extensions(
        path_folder: PathLike,
        extensions_to_keep: Tuple[str],
) -> List[str]:
    """List all the files with the provided extensions
    in the path_folder.

    Parameters
    ----------
    path_folder: Starting point for the tree listing.
    extensions_to_keep: Files with these extensions will have their
        hashes computed and tracked.

    Returns
    -------
    all_files: List of files with the correct extensions.
    """
    all_files = []
    for subdir, dirs, files in os.walk(path_folder):
        files.sort()
        for file in files:
            if file.lower().endswith(extensions_to_keep):
                all_files.append(os.path.join(subdir, file))
    return all_files


def create_list_hashes(
        path_folder: PathLike,
        extensions_to_keep: Tuple[str] = (".nii.gz", ".tsv", ".json"),
) -> Dict:
    """Computes a dictionary of files with their corresponding hashes.

    Parameters
    ----------
    path_folder: Starting point for the tree listing.
    extensions_to_keep: Files with these extensions will have their
        hashes computed and tracked.

    Returns
    -------
    dict
        Dictionary of the form {/path/to/file.extension: hash(file.extension)}
    """
    import hashlib

    def file_as_bytes(file):
        with file:
            return file.read()

    return {
        fname[len(str(path_folder)) :]: str(
            hashlib.md5(file_as_bytes(open(fname, "rb"))).digest()
        )
        for fname in list_files_with_extensions(path_folder, extensions_to_keep)
    }


def compare_folders_structures(
        path_folder: PathLike,
        list_hashes: PathLike,
):
    """Compares the structure of a folder against a reference.

    Parameters
    ----------
    path_folder: Starting point for the tree listing.
    list_hashes: Path to the pickled hash dictionary.
        The dictionary is assumed to be of the form
        {/path/to/file.extension: hash(file.extension)}
    """
    import pickle

    hashes_check = pickle.load(open(list_hashes, "rb"))
    hashes_new = create_list_hashes(path_folder)

    if hashes_check != hashes_new:
        error_message1 = ""
        error_message2 = ""
        for key in hashes_check:
            if key not in hashes_new:
                error_message1 += f"{key} not found !\n"
        for key in hashes_new:
            if key not in hashes_check:
                error_message2 += f"{key}'s creation was not expected !\n"
        raise ValueError(error_message1 + error_message2)
