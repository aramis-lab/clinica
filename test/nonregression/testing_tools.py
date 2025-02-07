# coding: utf8

import os
from enum import Enum
from os import PathLike
from pathlib import Path
from typing import Callable, Dict, Union

import numpy as np
import pandas as pd


def configure_paths(
    base_dir: Path,
    tmp_path: Path,
    name: str,
) -> tuple[Path, Path, Path]:
    """Configure paths for tests."""
    input_dir = base_dir / name / "in"
    ref_dir = base_dir / name / "ref"
    tmp_out_dir = tmp_path / name / "out"
    tmp_out_dir.mkdir(parents=True, exist_ok=False)

    return input_dir, tmp_out_dir, ref_dir


def likeliness_measure(
    file1: PathLike,
    file2: PathLike,
    threshold1: tuple,
    threshold2: tuple,
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
    data1 = nib.load(str(file1)).get_fdata(dtype="float32")
    data1[np.isnan(data1)] = 0

    data2 = nib.load(str(file2)).get_fdata(dtype="float32")
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
            [np.sum((metric_flattened > t)) / metric_flattened.size for t in thresholds]
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
    """Compare two NIfTI inputs using a similarity metric.

    Both NIfTI inputs are considered equal if the computed similarity metric
    is higher than the specified threshold.

    Parameters
    ----------
    file1: Path to first nifti input
    file2: Path to second nifti to compare
    threshold: Threshold value to be used in the comparison

    Returns
    -------
    bool
        True if file1 and file2 can be considered similar enough, i.e. the
        similarity metric is higher than the threshold.
    """
    from os import fspath

    import nibabel

    image1 = nibabel.load(fspath(file1)).get_fdata()
    image2 = nibabel.load(fspath(file2)).get_fdata()

    return similarity_measure_arrays(image1, image2, threshold)


def similarity_measure_large_nifti(
    file1: PathLike,
    file2: PathLike,
    threshold: float,
) -> bool:
    """Compare two NIfTI inputs using a similarity metric.

    Both NIfTI inputs are considered equal if the computed similarity metric
    is higher than the specified threshold.

    The difference with the function `similarity_measure` is that the two
    nifti images are compared slice by slice (taken in the last dimension).
    This allows to compare large 4D volumes that wouldn't fit entirely in
    memory for example.

    Parameters
    ----------
    file1 : Path
        The path to first nifti input.

    file2 : Path
        The path to second nifti to compare.

    threshold : float
        Threshold value to be used in the comparison.

    Returns
    -------
    bool
        True if file1 and file2 can be considered similar enough, i.e. the
        similarity metric is higher than the threshold."""
    from os import fspath

    import nibabel

    # Note that nibabel does lazy loading so image1 and 2 are only proxies here
    image1 = nibabel.load(fspath(file1))
    image2 = nibabel.load(fspath(file2))

    for volume in range(image1.shape[-1]):
        if not similarity_measure_arrays(
            np.asarray(image1.dataobj[..., volume]),
            np.asarray(image2.dataobj[..., volume]),
            threshold,
        ):
            return False
    return True


def similarity_measure_arrays(
    array1: np.ndarray,
    array2: np.ndarray,
    threshold: float,
) -> bool:
    """Returns True if structural similarity between two arrays is larger than threshold.

    See https://scikit-image.org/docs/stable/api/skimage.metrics.html#skimage.metrics.structural_similarity
    """
    import numpy as np
    from skimage.metrics import structural_similarity

    array1 = np.squeeze(array1)
    array2 = np.squeeze(array2)

    similarity = structural_similarity(
        array1,
        array2,
        gaussian_weights=True,
        sigma=1.5,
        use_sample_covariance=False,
        data_range=array1.max() - array1.min(),
    )

    return similarity > threshold


def compare_subject_session_tsv(
    sub_ses_tsv_1: PathLike,
    sub_ses_tsv_2: PathLike,
) -> bool:
    """Ensures that both subject_session TSV files are describing the same list of sessions.

    Parameters
    ----------
    sub_ses_tsv_1 : Path
        The path to first subject session TSV file.

    sub_ses_tsv_2: Path
        The path to second subject session TSV file.

    Returns
    -------
    bool
        True if sub_ses_tsv_1 and sub_ses_tsv_2 contains the same sessions.
    """
    # The operation is performed both sides because is_included(list1, list2) != is_included(list2, list1)
    return _is_included(sub_ses_tsv_1, sub_ses_tsv_2) & _is_included(
        sub_ses_tsv_2, sub_ses_tsv_1
    )


def _is_included(sub_ses_tsv_1: PathLike, sub_ses_tsv_2: PathLike) -> bool:
    df1 = pd.read_csv(sub_ses_tsv_1, sep="\t")
    df2 = pd.read_csv(sub_ses_tsv_2, sep="\t")

    if list(df1.columns) != list(df2.columns):
        return False
    subjects1 = list(df1.participant_id)
    sessions1 = list(df1.session_id)
    subjects2 = list(df2.participant_id)
    sessions2 = list(df2.session_id)

    if len(subjects1) != len(subjects2):
        return False
    for i in range(len(subjects1)):
        current_sub = subjects1[i]
        current_ses = sessions1[i]
        # Compute all the indices in the second list corresponding to the current subject
        idx_same_sub = [j for j in range(len(subjects2)) if subjects2[j] == current_sub]
        if len(idx_same_sub) == 0:  # Current subject not found in
            return False
        ses_same_sub = [sessions2[idx] for idx in idx_same_sub]
        if current_ses not in ses_same_sub:
            return False
    return True


def _sort_subject_field(subjects: list, fields: list) -> list:
    """Helper function for `same_missing_modality_tsv`.
    Returns a sorted list of fields. The list is sorted by corresponding
    subject_id and by field_id if the subject_ids are equal.
    """
    return [list(_) for _ in zip(*sorted(zip(subjects, fields)))][1]


def _extract_modality_from_tsv(file: PathLike) -> Dict:
    """Helper function for `same_missing_modality_tsv`.

    Extract data and create lists for a given TSV file.
    Subjects are sorted in alphabetical order. The same
    permutation of element is applied on each column.
    """
    df = pd.read_csv(file, sep="\t")
    subjects = list(df.participant_id)
    fields = {
        "pet_AV45": list(df["pet_trc-18FAV45"]),
        "pet_FDG": list(df["pet_trc-18FFDG"]),
        "t1w": list(df.t1w),
        "func_task_rest": list(df["func_task-rest"]),
    }
    subjects_sorted = sorted(subjects)
    sorted_fields = {k: _sort_subject_field(subjects, v) for k, v in fields.items()}
    sorted_fields["subjects_sorted"] = subjects_sorted
    return sorted_fields


def compare_missing_modality_tsv(file1: PathLike, file2: PathLike) -> bool:
    """Compare 2 TSV files generated by the iotool ComputeMissingModalities.

    Only fields participant_id, pet, t1w, func_task - rest are compared. Line order does not matter.

    Parameters
    ----------
    file1 : Path
        The path to first tsv.

    file2 : Path
        The path to second tsv.

    Returns
    -------
    bool
        True if file1 and file2 contains the same information
    """
    mods = [_extract_modality_from_tsv(f) for f in [file1, file2]]
    return mods[0] == mods[1]


def compare_folders(outdir: PathLike, refdir: PathLike, tmp_path: PathLike) -> bool:
    from filecmp import cmp
    from pathlib import PurePath

    file_out = PurePath(tmp_path) / "file_out.txt"
    file_ref = PurePath(tmp_path) / "file_ref.txt"
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


def tree(dir_: PathLike, file_out: PathLike):
    """Creates a file (file_out) with a visual tree representing the file
    hierarchy at a given directory

    .. note::
        Does not display empty directories.

    """
    from pathlib import Path

    file_content = ""

    for path in sorted(Path(dir_).rglob("*")):
        if path.is_dir() and not any(path.iterdir()):
            continue
        depth = len(path.relative_to(dir_).parts)
        spacer = "    " * depth
        file_content = file_content + f"{spacer}+ {path.name}\n"

    print(file_content)

    Path(file_out).write_text(file_content)


def clean_folder(path: PathLike, recreate: bool = True):
    """Clean folders under provided path.

    Parameters
    ----------
    path: Path of folder to clean.
    recreate: Whether to restore folder structure after cleaning.
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
    extensions_to_keep: tuple[str, ...],
) -> list[str]:
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
    extensions_to_keep: tuple[str, ...] = (".nii.gz", ".tsv", ".json"),
) -> Dict:
    """Computes a dictionary of files with their corresponding hashes.

    Parameters
    ----------
    path_folder : Path
        Starting point for the tree listing.

    extensions_to_keep : tuple of str
        Files with these extensions will have their hashes computed and tracked.

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


def compare_folders_with_hashes(
    path_folder: PathLike,
    list_hashes: PathLike,
):
    """Compares the files of a folder against a reference.

    Parameters
    ----------
    path_folder: Starting point for the tree listing.
    list_hashes: Path to the pickled hash dictionary.
        The dictionary is assumed to be of the form
        {/path/to/file.extension: hash(file.extension)}

    See Also
    --------
    compare_folders_structures
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
            elif hashes_check[key] != hashes_new[key]:
                error_message2 += f"{key} does not match the reference file !\n"
        raise ValueError(error_message1 + error_message2)


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

    See Also
    --------
    compare_folders_with_hashes
    """
    import pickle

    hashes_check = pickle.load(open(list_hashes, "rb"))
    hashes_new = create_list_hashes(path_folder)

    if set(hashes_check.keys()) != set(hashes_new.keys()):
        error_message1 = ""
        error_message2 = ""
        for key in hashes_check:
            if key not in hashes_new:
                error_message1 += f"{key} not found !\n"
        for key in hashes_new:
            if key not in hashes_check:
                error_message2 += f"{key}'s creation was not expected !\n"
        raise ValueError(error_message1 + error_message2)


class Level(str, Enum):
    PARTICIPANTS = "participants"
    SESSIONS = "sessions"
    SCANS = "scans"


def _load_participants_tsv(
    bids_dir: Path,
    _: Path,
) -> pd.DataFrame:
    return pd.read_csv(bids_dir / "participants.tsv", sep="\t").sort_values(
        by="participant_id", ignore_index=True
    )


def _load_sessions_tsv(bids_dir: Path, ref_tsv: Path) -> pd.DataFrame:
    return pd.read_csv(
        bids_dir / ref_tsv.parent.name / ref_tsv.name, sep="\t"
    ).sort_values(by="session_id", ignore_index=True)


def _load_scans_tsv(bids_dir: Path, ref_tsv: Path) -> pd.DataFrame:
    return pd.read_csv(
        bids_dir / ref_tsv.parent.parent.name / ref_tsv.parent.name / ref_tsv.name,
        sep="\t",
    ).sort_values(by="filename", ignore_index=True)


LoaderInterface = Callable[[Path, Path], pd.DataFrame]


def _loader_factory(level: Union[str, Level]) -> LoaderInterface:
    if (level := Level(level)) == Level.PARTICIPANTS:
        return _load_participants_tsv
    if level == Level.SESSIONS:
        return _load_sessions_tsv
    if level == Level.SCANS:
        return _load_scans_tsv
    raise (ValueError, f"TSV metadata file loader not implemented for level {level}.")


def _compare_frames(df1: pd.DataFrame, df2: pd.DataFrame, object: str):
    from pandas.testing import assert_frame_equal

    assert_frame_equal(df1, df2, check_like=True, obj=object)


def _iteratively_compare_frames(bids_ref: Path, bids_out: Path, level: Level):
    loader = _loader_factory(level)
    for tsv in bids_ref.rglob(f"*{level.value}.tsv"):
        _compare_frames(loader(bids_out, tsv), loader(bids_ref, tsv), tsv.name)


def compare_bids_tsv(bids_out: Path, bids_ref: Path):
    errors = []
    for level in Level:
        try:
            _iteratively_compare_frames(bids_ref, bids_out, level)
        except AssertionError as e:
            errors += [str(e).replace("\n\n", "\n")]
    if errors:
        raise AssertionError("\n\n".join(errors))
