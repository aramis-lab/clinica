import os.path
from os import PathLike
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple, Union

import nibabel as nib
import numpy as np
import pandas as pd

from clinica.pipelines.utils import FreeSurferAnnotationImage
from clinica.utils.image import HemiSphere
from clinica.utils.pet import SUVRReferenceRegion, Tracer

__all__ = [
    "get_output_dir",
    "remove_nan_from_image",
    "perform_gtmseg",
    "make_label_conversion",
    "run_apply_inverse_deformation_field",
    "run_apply_inverse_deformation_field_spm_standalone",
    "normalize_suvr",
    "reformat_surfname",
    "run_mris_expand",
    "run_mri_surf2surf",
    "run_mri_vol2surf",
    "compute_weighted_mean_surface",
    "project_onto_fsaverage",
    "get_mid_surface",
    "compute_average_pet_signal_based_on_annotations",
    "get_regexp_substitutions",
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

    from clinica.utils.filemanip import copy_file

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
    copy_file(gtmseg_filename, Path.cwd() / "gtmseg.mgz")

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
    """Method used on the segmentation from gtmsegmentation.

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


def run_apply_inverse_deformation_field(
    target_image: Path,
    deformation_field: Path,
    image: Path,
    matscript_folder: Path,
) -> Path:
    import sys

    from nipype.interfaces.matlab import MatlabCommand, get_matlab_command

    prefix = "subject_space_"
    MatlabCommand.set_default_matlab_cmd(get_matlab_command())
    matlab = MatlabCommand()
    if sys.platform.startswith("linux"):
        matlab.inputs.args = "-nosoftwareopengl"
    matlab.inputs.paths = str(matscript_folder)
    matlab.inputs.script = """
    applyInverseDeformationField('%s', '%s', '%s', './', '%s')
    """ % (
        target_image,
        deformation_field,
        image,
        prefix,
    )
    matlab.inputs.single_comp_thread = False
    matlab_log_file = Path.cwd() / "matlab_output.log"
    matlab.inputs.logfile = str(matlab_log_file)
    matlab.run()
    output_file = Path.cwd() / f"{prefix}{image.name}"
    if not output_file.exists():
        raise IOError(
            f"Something went wrong, please check {matlab_log_file.resolve()} for more information."
        )
    return output_file


def run_apply_inverse_deformation_field_spm_standalone(
    target_image: Path,
    deformation_field: Path,
    image: Path,
    matscript_folder: Path,
) -> Path:
    """Perform the same job as run_apply_inverse_deformation_field but with SPM standalone.

    We directly create a batch file that SPM standalone can run.
    This function does not check whether SPM standalone must be used.
    Previous check when building the pipeline ensures that all the
    env vars exists ($SPMSTANDALONE_HOME and $MCR_HOME).
    """
    import platform
    import subprocess

    prefix = "subject_space_"
    script_location = Path.cwd() / "m_script.m"
    with open(script_location, "w+") as script_file:
        script_file.write(
            _build_spm_batch_command(target_image, deformation_field, image, prefix)
        )
    if platform.system() == "Darwin":
        cmdline = f"cd $SPMSTANDALONE_HOME && ./run_spm12.sh $MCR_HOME batch {script_location}"
    elif platform.system() == "Linux":
        cmdline = f"$SPMSTANDALONE_HOME/run_spm12.sh $MCR_HOME batch {script_location}"
    else:
        raise SystemError("Only support Mac OS and Linux")
    subprocess_run = subprocess.run(
        cmdline,
        shell=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    if subprocess_run.returncode != 0:
        raise ValueError(
            "runApplyInverseDeformationField_SPM_standalone failed, returned non-zero code"
        )
    output_file = Path.cwd() / f"{prefix}{image.name}"
    if output_file.exists():
        raise IOError(
            "Something went wrong while trying to run runApplyInverseDeformationField_SPM_standalone. "
            "Output file not generated. Command launched :\n\t {cmdline}\n. We strongly recommend that "
            "you use the supported version of Matlab MCR recommended by the creators of SPM."
        )
    return output_file


def _build_spm_batch_command(
    target_image: Path, deformation_field: Path, image: Path, prefix: str
) -> str:
    """Write SPM batch command directly in a script that is readable by SPM standalone."""
    import os
    from os.path import abspath

    command = (
        "jobs{1}.spm.util.defs.comp{1}.inv.comp{1}.def = {'"
        + str(deformation_field)
        + "'};\n"
    )
    command += (
        "jobs{1}.spm.util.defs.comp{1}.inv.space = {'" + str(target_image) + "'};\n"
    )
    command += "jobs{1}.spm.util.defs.out{1}.pull.fnames = {'" + str(image) + "'};\n"
    command += (
        "jobs{1}.spm.util.defs.out{1}.pull.savedir.saveusr = {'"
        + abspath(os.getcwd())
        + "'};\n"
    )
    command += "jobs{1}.spm.util.defs.out{1}.pull.interp = 4;\n"
    command += "jobs{1}.spm.util.defs.out{1}.pull.mask = 1;\n"
    command += "jobs{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];\n"
    command += "jobs{1}.spm.util.defs.out{1}.pull.prefix = '" + prefix + "';\n"

    return command


def normalize_suvr(
    pet_image: Path, mask: Path, output_dir: Optional[Path] = None
) -> Path:
    """Get SUVR from pet image.

    Based on the segmentation performed by gtmsegmentation.
    The Standard Uptake Value ratio is computed by dividing the
    whole PET volume by the mean value observed in the pons.

    Parameters
    ----------
    pet_image : Path
        The path to the Nifti volume containing PET scan, realigned on up-sampled T1.

    mask : Path
        The path to the mask of the pons (18FFDG) or pons+cerebellum (18FAV45) already eroded.

    output_dir : Path, optional
        The directory in which to write the SUVR image.
        If not provided, it will be written in the current directory.

    Returns
    -------
    Path :
        The path to the SUVR normalized volume in the current directory.

    Raises
    ------
    ClinicaImageError :
        If the provided eroded mask contains only zero values.
    """
    import nibabel as nib

    from clinica.utils.exceptions import ClinicaImageError

    eroded_mask_nifti = nib.load(mask)
    eroded_mask = eroded_mask_nifti.get_fdata(dtype="float32")
    eroded_mask = eroded_mask > 0
    if (mask_size := np.sum(eroded_mask)) == 0:
        raise ClinicaImageError(
            f"The eroded mask located at {mask} contains only zero values. "
            "A problem likely occurred when moving the eroded mask from MNI to gtmsegspace."
        )
    # Load PET data (they must be in gtmsegspace, or same space as label file)
    pet_image_nifti = nib.load(pet_image)
    pet_data = pet_image_nifti.get_fdata(dtype="float32")
    # Mask unwanted values to determine mean uptake value
    pons_pet_activity = eroded_mask * pet_data
    mean_pons_pet_activity = np.sum(pons_pet_activity) / mask_size
    # Then normalize PET data by this mean activity
    suvr_image_nifti = nib.Nifti1Image(
        pet_data / mean_pons_pet_activity,
        pet_image_nifti.affine,
        header=pet_image_nifti.header,
    )
    suvr_image = (output_dir or Path.cwd()) / f"suvr_{pet_image.name}"
    nib.save(suvr_image_nifti, suvr_image)

    return suvr_image


def reformat_surfname(
    hemisphere: HemiSphere, left_surface: Path, right_surface: Path
) -> Path:
    if hemisphere == HemiSphere.LEFT:
        return left_surface
    if hemisphere == HemiSphere.RIGHT:
        return right_surface


def run_mris_expand(surface: Path, output_dir: Optional[Path] = None) -> List[Path]:
    """Make a subprocess call to the freesurfer mris_expand function.

    Expands the white input surface toward the pial, generating 7 surfaces at
    35%, 40%, 45%, 50%, 55%, 60%, 65% of thickness.

    Parameters
    ----------
    surface : Path
        The path to the input white surface.
        Must be named 'lh.white' or 'rh.white'.
        The folder containing the surface file must also have
        '?h.pial', '?.sphere', '?h.thickness' (freesurfer 'surf' folder).

    output_dir : Path, optional
        The path to the output folder in which to write the output files.
        If not provided, the files will be written in the current directory.

    Returns
    -------
    List of Path :
        List of path to the generated surfaces.

    Notes
    -----
    There is a bug in mris_expand : you are not allowed to write the surfaces
    elsewhere than in the surf folder.
    -N is a hidden parameter (not documented) that allows the user to specify
    the number of surface generated between source and final target surface.
    Here target is 65% of thickness, with 13 surfaces.
    Then we only keep the surfaces we are interested in.
    """
    from clinica.utils.filemanip import move_file
    from clinica.utils.stream import cprint

    output_dir = output_dir or Path.cwd()
    if not output_dir.exists():
        output_dir.mkdir(parents=True, exist_ok=True)
    out_file = surface.parent / f"{surface.name}_exp-"
    _run_mri_expand_as_subprocess(surface, out_file)
    # Move surface file to the output directory
    source_files = [
        out_file.parent / f"{out_file.name}0{n}"
        for n in ("07", "08", "09", "10", "11", "12", "13")
    ]
    destination_files = [output_dir / f.name for f in source_files]
    for src, dst in zip(source_files, destination_files):
        move_file(src, dst)
    # Remove useless surfaces (0%, 5%, 10%, 15%, 20%, 25% and 30% of thickness)
    for file_to_remove in (out_file.parent / f"{out_file.name}00{n}" for n in range(7)):
        file_to_remove.unlink(missing_ok=False)
        cprint(f"File {file_to_remove} has been removed.", lvl="debug")

    return destination_files


def _run_mri_expand_as_subprocess(surface: Path, out_file: Path):
    _run_command_as_subprocess(
        "mris_expand", _build_mri_expand_command(surface, out_file)
    )


def _run_command_as_subprocess(command_name: str, command: str):
    import subprocess

    from clinica.utils.exceptions import ClinicaSubprocessError
    from clinica.utils.stream import cprint

    cprint(
        f"Running {command_name} with the following command:\n\n{command}", lvl="debug"
    )
    subprocess_ = subprocess.run(
        command,
        shell=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    if (code := subprocess_.returncode) != 0:
        error_msg = (
            f"The subprocess '{command_name}' failed with a non-zero return code of {code}. "
            f"The following command was run:\n\n{command}"
        )
        cprint(error_msg, lvl="error")
        raise ClinicaSubprocessError(error_msg)


def _build_mri_expand_command(surface: Path, out_file: Path) -> str:
    import platform

    user_system = platform.system().lower()
    command = f"mris_expand -thickness -N 13 {surface} 0.65 {out_file}"
    # If system is MacOS, this export command must be run just
    # before the mri_vol2surf command to bypass MacOs security.
    if user_system.startswith("darwin"):
        command = _make_freesurfer_command_mac_compatible(command)
    return command


def _make_freesurfer_command_mac_compatible(command: str) -> str:
    return "export DYLD_LIBRARY_PATH=$FREESURFER_HOME/lib/gcc/lib && " + command


def run_mri_surf2surf(
    surface: Path,
    registration: Path,
    gtmsegfile: Path,
    subject_id: str,
    session_id: str,
    caps_dir: Path,
    is_longitudinal: bool,
    output_dir: Optional[Path] = None,
) -> Path:
    """Make a subprocess call to the freesurfer surf2surf function.

    Here the aim is to convert a input surface (which is the native space of the subject),
    into the same surface but in the gtmseg space (space of the volume generated by
    the gtmsegmentation).

    Parameters
    ----------
    surface : Path
        The path to the surface file that needs to be converted.

    registration : Path
        The path to a registration file that represents the transformation
        needed to go from the native space to the gtmsegspace
        (see https://surfer.nmr.mgh.harvard.edu/fswiki/FsAnat-to-NativeAnat
        for more details).

    gtmsegfile : Path
        The path to the gtm segmentation file.

    subject_id : str
        The subject_id (something like sub-ADNI002S4213).

    session_id : str
        The session id ( something like : ses-M012).

    caps_dir : Path
        The path to the CAPS directory.

    is_longitudinal : bool
        Whether the function should handle longitudinal files or not.

    output_dir : Path, optional
        The path to the output folder in which to write the output files.
        If not provided, the files will be written in the current directory.

    Returns
    -------
    Path :
        The path to the converted surface in current directory.
    """
    import os

    from clinica.utils.filemanip import copy_file

    # set subjects_dir env. variable for mri_surf2surf to work properly
    subjects_directory_backup = os.path.expandvars("$SUBJECTS_DIR")
    subjects_directory, freesurfer_id = (
        _get_new_subjects_directory_longitudinal
        if is_longitudinal
        else _get_new_subjects_directory
    )(caps_dir, subject_id, session_id)
    os.environ["SUBJECTS_DIR"] = str(subjects_directory)
    copy_file(surface, subjects_directory / freesurfer_id / "surf")
    output_file = (output_dir or Path.cwd()) / f"{surface.name}_gtmsegspace"
    _run_mri_surf2surf_as_subprocess(
        surface, registration, gtmsegfile, freesurfer_id, output_file
    )
    (subjects_directory / freesurfer_id / "surf" / surface.name).unlink(
        missing_ok=False
    )
    # put back original subjects_dir env
    os.environ["SUBJECTS_DIR"] = subjects_directory_backup

    return output_file


def _get_new_subjects_directory_longitudinal(
    caps_dir: Path,
    subject_id: str,
    session_id: str,
) -> Tuple[Path, str]:
    """Extract SUBJECT_DIR.

    Extract path to FreeSurfer segmentation in CAPS folder and FreeSurfer ID
    (e.g. sub-CLNC01_ses-M000.long.sub-CLNC01_long-M000M018 or sub-CLNC01_ses-M000).
    """
    from clinica.utils.exceptions import ClinicaCAPSError
    from clinica.utils.stream import cprint

    root = caps_dir / "subjects" / subject_id / session_id / "t1"
    long_folds = [f.name for f in root.iterdir() if f.name.startswith("long-")]
    if len(long_folds) > 1:
        error_msg = f"Folder {root} contains {len(long_folds)} folders labeled long-*. Only 1 can exist."
        cprint(error_msg, lvl="error")
        raise ClinicaCAPSError(error_msg)
    if len(long_folds) == 0:
        error_msg = (
            f"Folder {root} does not contains a folder labeled long-*. "
            "Have you run t1-freesurfer-longitudinal?"
        )
        cprint(error_msg, lvl="error")
        raise ClinicaCAPSError(error_msg)
    root_env = root / long_folds[0] / "freesurfer_longitudinal"
    sub_id_cmd = f"{subject_id}_{session_id}.long.{subject_id}_{long_folds[0]}"
    return root_env, sub_id_cmd


def _get_new_subjects_directory(
    caps_dir: Path, subject_id: str, session_id: str
) -> Tuple[Path, str]:
    return (
        caps_dir
        / "subjects"
        / subject_id
        / session_id
        / "t1"
        / "freesurfer_cross_sectional",
        f"{subject_id}_{session_id}",
    )


def _build_mri_surf2surf_command(
    surface: Path,
    registration: Path,
    gtmsegfile: Path,
    freesurfer_id: str,
    output_file: Path,
) -> str:
    import platform

    from clinica.utils.image import HemiSphere

    user_system = platform.system().lower()
    hemisphere = HemiSphere(surface.name[0:2])
    surface_name = surface.name[3:]
    command = (
        f"mri_surf2surf --reg {registration} {gtmsegfile} --sval-xyz {surface_name} "
        f"--hemi {hemisphere.value} --tval-xyz {gtmsegfile} --tval {output_file} --s {freesurfer_id}"
    )
    # If system is MacOS, this export command must be run just
    # before the mri_vol2surf command to bypass MacOs security.
    if user_system.startswith("darwin"):
        command = _make_freesurfer_command_mac_compatible(command)
    return command


def _run_mri_surf2surf_as_subprocess(
    surface: Path,
    registration: Path,
    gtmsegfile: Path,
    freesurfer_id: str,
    output_file: Path,
):
    _run_command_as_subprocess(
        "mri_surf2surf",
        _build_mri_surf2surf_command(
            surface, registration, gtmsegfile, freesurfer_id, output_file
        ),
    )


def run_mri_vol2surf(
    pet_volume: Path,
    surface: Path,
    subject_id: str,
    session_id: str,
    caps_dir: Path,
    gtmsegfile: Path,
    is_longitudinal: bool,
    output_dir: Optional[Path] = None,
) -> Path:
    """Make a subprocess call to the freesurfer vol2surf function.

    Projects the volume into the surface : the value at each vertex is
    given by the value of the voxel it intersects

    Parameters
    ----------
    pet_volume : Path
        The path to PET volume (in gtmseg space) that needs to be mapped into surface.

    surface : Path
        The path to surface file.

    subject_id : str
        The subject_id (something like sub-ADNI002S4213).

    session_id : str
        The session id ( something like : ses-M012).

    caps_dir : Path
        The path to the CAPS directory.

    gtmsegfile : Path
        The path to the gtm segmentation file (provides information on space, labels are not used).

    is_longitudinal : bool
        Whether the function should handle longitudinal files or not.

    output_dir : Path, optional
        The path to the output folder in which to write the output files.
        If not provided, the files will be written in the current directory.

    Returns
    -------
    Path :
        The path to the data projected onto the surface.
    """
    import os

    from clinica.utils.filemanip import copy_file
    from clinica.utils.image import HemiSphere

    subjects_directory_backup = os.path.expandvars("$SUBJECTS_DIR")
    subjects_directory, freesurfer_id = (
        _get_new_subjects_directory_longitudinal
        if is_longitudinal
        else _get_new_subjects_directory
    )(caps_dir, subject_id, session_id)
    os.environ["SUBJECTS_DIR"] = str(subjects_directory)
    copy_file(surface, subjects_directory / freesurfer_id / "surf")
    gtmsegfile_copy = subjects_directory / freesurfer_id / "mri" / "gtmseg.mgz"
    if not gtmsegfile_copy.exists():
        copy_file(gtmsegfile, gtmsegfile_copy)
    hemisphere = HemiSphere(surface.name[0:2])
    output_file = (
        output_dir or Path.cwd()
    ) / f"{hemisphere.value}.projection_{surface.name}.mgh"
    _run_mri_vol2surf_as_subprocess(pet_volume, surface, freesurfer_id, output_file)
    (subjects_directory / freesurfer_id / "surf" / surface.name).unlink(
        missing_ok=False
    )
    # TODO careful here...
    # Removing gtmseg.mgz may lead to problems as other vol2surf are using it
    gtmsegfile_copy.unlink(missing_ok=False)
    # put back original subjects_dir env
    os.environ["SUBJECTS_DIR"] = subjects_directory_backup

    return output_file


def _run_mri_vol2surf_as_subprocess(
    pet_volume: Path,
    surface: Path,
    freesurfer_id: str,
    output_file: Path,
):
    _run_command_as_subprocess(
        "mri_vol2surf",
        _build_mri_vol2surf_command(pet_volume, surface, freesurfer_id, output_file),
    )


def _build_mri_vol2surf_command(
    pet_volume: Path,
    surface: Path,
    freesurfer_id: str,
    output_file: Path,
) -> str:
    import platform

    from clinica.utils.image import HemiSphere

    user_system = platform.system().lower()
    hemisphere = HemiSphere(surface.name[0:2])
    surface_name = surface.name[3:]
    command = (
        f"mri_vol2surf --mov {pet_volume} --o {output_file} --surf {surface_name} --hemi {hemisphere.value} "
        f"--regheader {freesurfer_id} --ref gtmseg.mgz --interp nearest"
    )
    if user_system.startswith("darwin"):
        command = _make_freesurfer_command_mac_compatible(command)

    return command


def compute_weighted_mean_surface(
    surfaces: Sequence[Path], output_dir: Optional[Path] = None
) -> Path:
    """Compute a weighted average at each node of the surface.

    The weight are defined by a normal distribution (centered on the mid surface).

    Parameters
    ----------
    surfaces : Sequence of Path
        The paths to the data projected on the 7 surfaces (35 to 65 % of thickness) at each nodes.

    output_dir : Path, optional
        The path to the output folder in which to write the output files.
        If not provided, the files will be written in the current directory.

    Returns
    -------
    Path :
        The path to the data averaged.
    """
    import nibabel as nib

    _assert_seven_surfaces(surfaces)
    hemisphere = HemiSphere(surfaces[0].name[0:2])
    out_surface = (
        output_dir or Path.cwd()
    ) / f"{hemisphere.value}.averaged_projection_on_cortical_surface.mgh"
    nib.save(_build_weighted_mean_surface_image(surfaces), out_surface)

    return out_surface


def _assert_seven_surfaces(surfaces: Sequence[Path]):
    if (n_surfaces := len(surfaces)) != 7:
        raise ValueError(
            "There should be 7 surfaces at this point of the pipeline. "
            f"However 'compute_weighted_mean_surface' received {n_surfaces} surfaces. "
            "Something probably went wrong in prior steps of the pipeline."
        )


def _build_weighted_mean_surface_image(surfaces: Sequence[Path]) -> nib.MGHImage:
    import nibabel as nib
    import numpy as np

    # sample only to get dimension
    sample = nib.load(surfaces[0])
    data_normalized = np.zeros(sample.header.get_data_shape())
    for surface, coefficient in zip(
        surfaces, _get_coefficient_for_normal_repartition()
    ):
        current_surf = nib.load(surface)
        data_normalized += current_surf.get_fdata(dtype="float32") * coefficient
    # data_normalized = np.atleast_3d(data_normalized)
    return nib.MGHImage(data_normalized, affine=sample.affine, header=sample.header)


def _get_coefficient_for_normal_repartition() -> (
    Tuple[float, float, float, float, float, float, float]
):
    """TODO: Find out where do these numbers come from.."""
    return 0.1034, 0.1399, 0.1677, 0.1782, 0.1677, 0.1399, 0.1034


def project_onto_fsaverage(
    projection: Path,
    subject_id: str,
    session_id: str,
    caps_dir: Path,
    fwhm: int,
    is_longitudinal: bool,
    output_dir: Optional[Path] = None,
) -> Path:
    """Project data onto an averaged subject called fsaverage.

    This subject is available in the $SUBJECTS_DIR folder.

    Notes
    -----
    fsaverage and the subject must be in the subject_dir, so a copy of fsaverage is performed if necessary.

    Parameters
    ----------
    projection : Path
        The path to the projected data onto native subject surface.

    subject_id : str
        The subject id (something like sub-ADNI002S4213).

    session_id : str
        The session id ( something like : ses-M012).

    caps_dir : Path
        The path to the CAPS directory.

    fwhm : int
        FWHM of the Gaussian filter used for smoothing on fsaverage surface (not volume !)

    is_longitudinal : bool
        Longitudinal pipeline or not.

    output_dir : Path, optional
        The path to the output folder in which to write the output files.
        If not provided, the files will be written in the current directory.

    Returns
    -------
    Path :
        The path to the data averaged.
    """
    import os
    import shutil

    subjects_directory_backup = Path(os.path.expandvars("$SUBJECTS_DIR"))
    subjects_directory, freesurfer_id = (
        _get_new_subjects_directory_longitudinal
        if is_longitudinal
        else _get_new_subjects_directory
    )(caps_dir, subject_id, session_id)
    os.environ["SUBJECTS_DIR"] = str(subjects_directory)
    # copy fsaverage folder next to : subject_id + '_' + session_id
    # for the mris_preproc command to properly find src and target
    fsaverage_has_been_copied = False
    if not (subjects_directory / "fsaverage").exists():
        shutil.copytree(
            str(subjects_directory_backup / "fsaverage"),
            str(subjects_directory / "fsaverage"),
        )
        fsaverage_has_been_copied = True
    # also copy the mgh file in the surf folder (needed by MRISPreproc
    projection_in_surf_folder = (
        subjects_directory / freesurfer_id / "surf" / projection.name
    )
    if not projection_in_surf_folder.exists():
        shutil.copy(str(projection), str(projection_in_surf_folder))
    out_fsaverage = (
        output_dir or Path.cwd()
    ) / f"fsaverage_fwhm-{fwhm}_{projection.name}"
    _run_mris_preproc_as_standalone_nipype_node(
        projection, freesurfer_id, fwhm, out_fsaverage
    )
    # remove projection file from surf folder
    projection_in_surf_folder.unlink(missing_ok=False)
    # remove fsaverage if it has been copied
    if fsaverage_has_been_copied:
        shutil.rmtree(subjects_directory / "fsaverage")
    # put back original subjects_dir env
    os.environ["SUBJECTS_DIR"] = str(subjects_directory_backup)
    return out_fsaverage


def _run_mris_preproc_as_standalone_nipype_node(
    projection: Path,
    freesurfer_id: str,
    fwhm: float,
    output_file: Path,
):
    from nipype.interfaces.freesurfer import MRISPreproc

    projection_node = MRISPreproc()
    projection_node.inputs.target = "fsaverage"
    projection_node.inputs.subjects = [freesurfer_id]
    projection_node.inputs.fwhm = fwhm
    projection_node.inputs.hemi = HemiSphere(projection.name[0:2]).value
    projection_node.inputs.surf_measure = projection.name[3:]
    projection_node.inputs.out_file = str(output_file)
    projection_node.run()


def get_mid_surface(surfaces: Sequence[Path]) -> Path:
    """Returns the mid-surface when dealing with the 7 different surfaces.

    Parameters
    ----------
    surfaces : Sequence of Path
        The 7 different surfaces generated by mris_expand.

    Returns
    -------
    Path :
        The path to the mid-surface.
    """
    _assert_seven_surfaces(surfaces)
    return surfaces[3]


def compute_average_pet_signal_based_on_annotations(
    pet_projections: Tuple[Path, Path],
    atlas_files: Dict[str, FreeSurferAnnotationImage],
    output_dir: Optional[Path] = None,
) -> List[Path]:
    """Computes the average of PET signal based on annot files from Freesurfer.

    Those files describe the brain according to known atlases.

    Parameters
    ----------
    pet_projections : tuple of two Path
        The paths to the PET projection (must be a MGH file) [left_hemisphere, right_hemisphere].

    atlas_files : dict[str, FreeSurferAnnotationImage]
        Dictionary containing path to lh and rh annotation files for any number of atlases.

    output_dir : Path, optional
        The path to the output folder in which to write the output files.
        If not provided, the files will be written in the current directory.

    Returns
    -------
    Path :
        The path to the tsv containing average PET values.

    Raises
    ------
    ValueError :
        If not exactly two files were provided through the argument 'pet_projections'.
    """
    import nibabel as nib
    import numpy as np
    import pandas as pd

    from clinica.pipelines.utils import FreeSurferAnnotation
    from clinica.utils.stream import log_and_raise

    if len(pet_projections) != 2:
        msg = (
            "The compute_average_pet_signal_based_on_annotations function requires two files "
            "for the argument 'pet_projections', one for the left hemisphere, one for the right. "
            f"The following {len(pet_projections)} were received:\n"
            + "\n".join([str(_) for _ in pet_projections])
        )
        log_and_raise(msg, ValueError)
    pet_mgh = {
        HemiSphere.LEFT: np.squeeze(
            nib.load(pet_projections[0]).get_fdata(dtype="float32")
        ),
        HemiSphere.RIGHT: np.squeeze(
            nib.load(pet_projections[1]).get_fdata(dtype="float32")
        ),
    }
    filename_tsv = []
    for atlas_name, annotation_image in atlas_files.items():
        annotation = FreeSurferAnnotation.from_annotation_image(annotation_image)
        annotation.replace_minus_one_annotation_with_zero()
        average_region = []
        for region_id, region_name in enumerate(annotation.region_names):
            for hemisphere in (HemiSphere.LEFT, HemiSphere.RIGHT):
                mask = annotation.get_annotation(hemisphere) == region_id
                mask = np.uint(mask)
                masked_data = mask * pet_mgh[hemisphere]
                average_region.append(
                    np.nan if np.sum(mask) == 0 else np.sum(masked_data) / np.sum(mask)
                )
        final_tsv = pd.DataFrame(
            {
                "index": range(len(average_region)),
                "label_name": annotation.get_lateralized_region_names(left_first=True),
                "mean_scalar": average_region,
            }
        )
        filename_atlas_tsv = (output_dir or Path.cwd()) / f"{atlas_name}.tsv"
        filename_tsv.append(filename_atlas_tsv)
        final_tsv.to_csv(
            filename_atlas_tsv,
            sep="\t",
            index=False,
            columns=["index", "label_name", "mean_scalar"],
        )
    return filename_tsv


def get_regexp_substitutions(
    pet_tracer: Tracer,
    region: SUVRReferenceRegion,
    is_longitudinal: bool,
) -> List[Tuple[str, str]]:
    return [
        _get_mid_surface_substitutions(is_longitudinal=is_longitudinal),
        _get_projection_in_native_space_substitutions(
            pet_tracer, region, is_longitudinal=is_longitudinal
        ),
        _get_projection_in_fsaverage_substitution(
            pet_tracer, region, is_longitudinal=is_longitudinal
        ),
        _get_tsv_file_for_atlas(
            pet_tracer, region, "destrieux", is_longitudinal=is_longitudinal
        ),
        _get_tsv_file_for_atlas(
            pet_tracer, region, "desikan", is_longitudinal=is_longitudinal
        ),
    ]


def _get_mid_surface_substitutions(is_longitudinal: bool) -> Tuple[str, str]:
    if is_longitudinal:
        return (
            r"(.*(sub-.*)\/(ses-.*)\/pet\/(long-.*)\/surface_longitudinal)\/midsurface\/.*_hemi_([a-z]+)(.*)$",
            r"\1/\2_\3_\4_hemi-\5_midcorticalsurface",
        )
    return (
        r"(.*(sub-.*)\/(ses-.*)\/pet\/surface)\/midsurface\/.*_hemi_([a-z]+)(.*)$",
        r"\1/\2_\3_hemi-\4_midcorticalsurface",
    )


def _get_projection_in_native_space_substitutions(
    pet_tracer: Tracer,
    region: SUVRReferenceRegion,
    is_longitudinal: bool,
) -> Tuple[str, str]:
    if is_longitudinal:
        return (
            r"(.*(sub-.*)\/(ses-.*)\/pet\/(long-.*)\/surface_longitudinal)\/projection_native\/.*_hemi_([a-z]+).*",
            rf"\1/\2_\3_\4_trc-{pet_tracer.value}_pet_space-native_suvr-{region.value}_pvc-iy_hemi-\5_projection.mgh",
        )
    return (
        r"(.*(sub-.*)\/(ses-.*)\/pet\/surface)\/projection_native\/.*_hemi_([a-z]+).*",
        rf"\1/\2_\3_trc-{pet_tracer.value}_pet_space-native_suvr-{region.value}_pvc-iy_hemi-\4_projection.mgh",
    )


def _get_projection_in_fsaverage_substitution(
    pet_tracer: Tracer,
    region: SUVRReferenceRegion,
    is_longitudinal: bool,
) -> Tuple[str, str]:
    if is_longitudinal:
        return (
            (
                r"(.*(sub-.*)\/(ses-.*)\/pet\/(long-.*)\/surface_longitudinal)\/"
                r"projection_fsaverage\/.*_hemi_([a-z]+).*_fwhm_([0-9]+).*"
            ),
            (
                rf"\1/\2_\3_\4_trc-{pet_tracer.value}_pet_space-fsaverage_"
                rf"suvr-{region.value}_pvc-iy_hemi-\5_fwhm-\6_projection.mgh"
            ),
        )
    return (
        r"(.*(sub-.*)\/(ses-.*)\/pet\/surface)\/projection_fsaverage\/.*_hemi_([a-z]+).*_fwhm_([0-9]+).*",
        (
            rf"\1/\2_\3_trc-{pet_tracer.value}_pet_space-fsaverage_"
            rf"suvr-{region.value}_pvc-iy_hemi-\4_fwhm-\5_projection.mgh"
        ),
    )


def _get_tsv_file_for_atlas(
    pet_tracer: Tracer,
    region: SUVRReferenceRegion,
    atlas: str,
    is_longitudinal: bool,
) -> Tuple[str, str]:
    if is_longitudinal:
        return (
            rf"(.*(sub-.*)\/(ses-.*)\/pet\/(long-.*)\/surface_longitudinal)\/{atlas}_tsv\/{atlas}.tsv",
            (
                rf"\1/atlas_statistics/\2_\3_\4_trc-{pet_tracer.value}_pet_"
                rf"space-{atlas}_pvc-iy_suvr-{region.value}_statistics.tsv"
            ),
        )
    return (
        rf"(.*(sub-.*)\/(ses-.*)\/pet\/surface)\/{atlas}_tsv\/{atlas}.tsv",
        (
            rf"\1/atlas_statistics/\2_\3_trc-{pet_tracer.value}_pet_"
            rf"space-{atlas}_pvc-iy_suvr-{region.value}_statistics.tsv"
        ),
    )
