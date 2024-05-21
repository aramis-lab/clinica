import os.path
from os import PathLike
from pathlib import Path
from typing import List, Tuple, Union

import nibabel as nib
import numpy as np
import pandas as pd

__all__ = [
    "get_output_dir",
    "remove_nan_from_image",
    "perform_gtmseg",
]


def get_output_dir(
    is_longitudinal: bool,
    caps_dir: Union[str, PathLike],
    subject_id: str,
    session_id: str,
) -> Path:
    caps_dir = Path(caps_dir)
    root = caps_dir / "subjects" / subject_id / session_id
    if is_longitudinal:
        input_folder = root / "t1"
        return (
            root
            / "pet"
            / _get_longitudinal_folder_name(input_folder)
            / "surface_longitudinal"
        )
    return root / "pet" / "surface"


def _get_longitudinal_folder_name(input_folder: Path) -> str:
    from clinica.utils.exceptions import ClinicaCAPSError

    longitudinal_folders = [
        f.name for f in input_folder.iterdir() if f.name.startswith("long-")
    ]
    if len(longitudinal_folders) > 1:
        raise ClinicaCAPSError(
            f"[Error] Folder {input_folder} contains {len(longitudinal_folders)} "
            "folders labeled long-*. Only 1 can exist"
        )
    if len(longitudinal_folders) == 0:
        raise ClinicaCAPSError(
            f"[Error] Folder {input_folder} does not contains a folder labeled long-*. "
            "Have you run t1-freesurfer-longitudinal?"
        )
    return longitudinal_folders[0]


def remove_nan_from_image(image_path: PathLike) -> Path:
    """Remove NaN values from the provided nifti image.

    This is needed after a registration performed by 'spmregister' : instead
    of filling space with 0, nan are used to extend the PET space.
    We propose to replace them with 0s.

    Parameters
    ----------
    image_path : PathLike
        The path to the Nifti volume where NaNs need to be replaced by zeros.

    Returns
    -------
    output_image_path : Path
        The path to the volume in Nifti that does not contain any NaNs.
    """
    import nibabel as nib
    import numpy as np

    from clinica.utils.filemanip import get_filename_no_ext

    image = nib.load(image_path)
    data = np.nan_to_num(image.get_fdata(dtype="float32"))
    output_image = nib.Nifti1Image(data, image.affine, header=image.header)
    output_image_path = Path.cwd() / f"no_nan_{get_filename_no_ext(image_path)}.nii.gz"
    nib.save(output_image, output_image_path)

    return output_image_path


def perform_gtmseg(
    caps_dir: Path, subject_id: str, session_id: str, is_longitudinal: bool
) -> Path:
    """Perform Freesurfer gtmseg.

    'gtmseg' is a freesurfer command used to perform a segmentation
    used in some partial volume correction methods.

    Parameters
    ----------
    caps_dir : Path
        CAPS directory.

    subject_id : str
        The subject ID. Example: 'sub-ADNI002S4213'.

    session_id : str
        The session ID. Example: 'ses-M012'.

    is_longitudinal : bool
        If longitudinal processing, subjects_dir must be put elsewhere

    Returns
    -------
    Path :
        Path to the segmentation volume : a volume where each voxel
        has a label (ranging [0 2035] ), see Freesurfer lookup table to see the
        labels with their corresponding names.

    Warnings
    --------
    This method changes the environment variable $SUBJECTS_DIR (but put
    the original one back after execution). This has not been intensely
    tested whether it can lead to some problems : (for instance if 2
    subjects are running in parallel)
    """
    import os
    import shutil

    # Old subject_dir is saved for later
    subjects_dir_backup = _expand_environment_variable_into_path("$SUBJECTS_DIR")
    root_env, freesurfer_id = _get_new_subjects_dir(
        is_longitudinal, caps_dir, subject_id, session_id
    )
    # Set the new subject dir for the function to work properly
    os.environ["SUBJECTS_DIR"] = str(root_env)

    freesurfer_mri_folder = (
        _expand_environment_variable_into_path("$SUBJECTS_DIR") / freesurfer_id / "mri"
    )
    gtmseg_filename = freesurfer_mri_folder / "gtmseg.mgz"
    if not gtmseg_filename.exists():
        _run_gtmseg(freesurfer_id)
    # We specify the out file to be in the current directory of execution
    # (easy for us to look at it afterward in the working directory).
    # We copy then the file.
    shutil.copy(gtmseg_filename, Path("gtmseg.mgz").resolve())

    # Remove bunch of files created during segmentation in caps dir and not needed
    for filename in ("gtmseg.ctab", "gtmseg.lta"):
        if (freesurfer_mri_folder / filename).exists():
            (freesurfer_mri_folder / filename).unlink()

    # Set back the SUBJECT_DIR environment variable of the user
    os.environ["SUBJECTS_DIR"] = str(subjects_dir_backup)

    return gtmseg_filename


def _expand_environment_variable_into_path(variable_name: str) -> Path:
    return Path(os.path.expandvars(variable_name))


def _get_new_subjects_dir(
    is_longitudinal: bool, caps_dir: Path, subject_id: str, session_id: str
) -> Tuple[Path, str]:
    """Extract SUBJECT_DIR.

    Extract path to FreeSurfer segmentation in CAPS folder and FreeSurfer ID
    (e.g. sub-CLNC01_ses-M000.long.sub-CLNC01_long-M000M018 or sub-CLNC01_ses-M000).
    """
    root = caps_dir / "subjects" / subject_id / session_id / "t1"
    if is_longitudinal:
        longitudinal_folder_name = _get_longitudinal_folder_name(root)
        return (
            root / longitudinal_folder_name / "freesurfer_longitudinal",
            f"{subject_id}_{session_id}.long.{subject_id}_{longitudinal_folder_name}",
        )
    return root / "freesurfer_cross_sectional", f"{subject_id}_{session_id}"


def _run_gtmseg(freesurfer_id: str):
    """Run the gtmseg command with provided freesurfer ID.

    This function creates a standalone node based on Command Line Interface.
    We simply put the command line we would run on a console.
    """
    import nipype.pipeline.engine as pe
    from nipype.interfaces.base import CommandLine

    segmentation = pe.Node(
        interface=CommandLine(
            f"gtmseg --s {freesurfer_id} --no-seg-stats --xcerseg",
            terminal_output="stream",
        ),
        name="gtmseg",
    )
    segmentation.run()


def make_label_conversion(gtmseg_file: Path, csv_file: Path) -> List[Path]:
    """

    'make_label_conversion' is a method used on the segmentation from gtmsegmentation.
    The purpose is to reduce the number of label.
    The gathering of labels is specified in a separate file.

    Parameters
    ----------
    gtmseg_file : Path
        The path to the Nifti volume containing the gtmseg segmentation.

    csv_file : Path
        The path to .csv file that contains 3 columns : REGION SOURCE DST.
        Separator is , (coma).

    Returns
    -------
    list of Path :
        List of path to the converted volumes according to the .csv file.
        Each volume is a mask representing an area.
    """
    label_image = nib.load(gtmseg_file)
    label_image.header.set_data_dtype("int8")
    label_data = label_image.get_fdata(dtype="float32")
    initial_labels = np.unique(label_data).astype("int16")
    control_volume = np.zeros(label_data.shape)
    src_val, dst_val, reg = _get_source_destination_region(csv_file)
    _verify_all_initial_labels_have_matching_transformation(initial_labels, src_val)

    # Instantiation of final volume, with same dtype as original volume
    new_label_data = np.zeros(label_data.shape, dtype=label_data.dtype)
    # Computing the transformation
    for i, src in enumerate(src_val):
        new_label_data[label_data == src] = dst_val[i]
    # Get unique list of new label
    new_labels = np.unique(new_label_data)
    new_labels = new_labels.astype("int")
    list_of_regions = []

    for label_value in new_labels:
        region_volume = np.zeros(label_data.shape, dtype="uint8")
        region_volume[new_label_data == label_value] = 1
        region_image = nib.Nifti1Image(
            region_volume, label_image.affine, header=label_image.header
        )
        output_path = Path(f"{label_value}.nii.gz").resolve()
        nib.save(region_image, output_path)
        list_of_regions.append(output_path)
        control_volume = control_volume + region_volume
    _check_sum(control_volume)

    return list_of_regions


def _read_csv(csv_file: Path) -> pd.DataFrame:
    expected_columns = ["REGION", "SOURCE", "DST"]
    if not csv_file.is_file():
        raise IOError(f"The provided CSV file {csv_file} does not exist.")
    df = pd.read_csv(csv_file, sep=",")
    if list(df.columns.values) != expected_columns:
        raise Exception(
            f"CSV file {csv_file} is not in the correct format. "
            f"Columns should be: {expected_columns}."
        )
    return df


def _verify_all_initial_labels_have_matching_transformation(
    initial_labels: np.ndarray, src_val: np.ndarray
):
    """Check that each label of original volume (initial_labels)
    has a matching transformation in the csv file.
    """
    for label in initial_labels:
        index = np.argwhere(src_val == label)
        # Size 0 means no occurrence found
        if index.size == 0:
            raise ValueError(
                f"Could not find label {label} on conversion table. "
                "Add it manually in CSV file to correct error."
            )


def _get_source_destination_region(
    csv_file: Path,
) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    df = _read_csv(csv_file)

    return (
        np.asanyarray(list(df.SOURCE)).astype("int"),
        np.asarray(list(df.DST)).astype("int"),
        list(df.REGION),
    )


def _check_sum(control_image_data: np.ndarray):
    """The sum of a voxel location across the fourth dimension should be 1."""
    sum_voxel_mean = float(sum(sum(sum(control_image_data)))) / control_image_data.size
    if not _are_almost_equal(1.0, sum_voxel_mean):
        raise ValueError(
            "Problem during parcellation: the mean sum of a voxel across "
            f"4th dimension is {sum_voxel_mean} instead of 1.0"
        )


def _are_almost_equal(a: float, b: float, rel_tol=1e-9, abs_tol=0.0) -> bool:
    """Measure equality between to floating or double numbers, using 2 thresholds : a
    relative tolerance, and an absolute tolerance
    """
    return abs(a - b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)
