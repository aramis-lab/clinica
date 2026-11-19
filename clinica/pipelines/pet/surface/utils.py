import os
import platform
import shutil
from pathlib import Path
from typing import Sequence

import nibabel as nib
import numpy as np
import pandas as pd

from clinica.pipelines.utils import FreeSurferAnnotation
from clinica.utils.filemanip import copy_file, move_file
from clinica.utils.image import HemiSphere
from clinica.utils.pet import SUVRReferenceRegion, Tracer
from clinica.utils.stream import cprint
from clinica.utils.third_party_execution import run_command_as_subprocess

__all__ = [
    "perform_gtmseg",
    "convert_labels",
    "run_mri_vol2surf",
    "compute_weighted_mean_surface",
    "project_onto_fsaverage",
    "get_mid_surface",
    "reformat_surfname",
    "compute_average_pet_signal_based_on_annotations",
    "merge_nifti_volumes",
    "run_apply_inverse_deformation_field_SPM_standalone",
    "run_mri_surf2surf",
    "normalize_suvr",
    "run_mris_expand",
    "remove_nan_from_image",
    "get_regexp_substitutions",
    "get_output_dir",
]


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


def get_output_dir(
    is_longitudinal: bool, caps_dir: Path, subject_id: str, session_id: str
) -> Path:
    root = caps_dir / "subjects" / subject_id / session_id
    if is_longitudinal:
        return (
            root
            / "pet"
            / _get_longitudinal_folder_name(root / "t1")
            / "surface_longitudinal"
        )

    return root / "pet" / "surface"


def _get_new_subjects_dir(
    is_longitudinal: bool,
    caps_dir: Path,
    subject_id: str,
    session_id: str,
) -> tuple[Path, str]:
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


def _run_gtmseg(freesurfer_id: str, subjects_directory: Path):
    """Run the gtmseg command with provided freesurfer ID.
    This function creates a standalone node based on Command Line Interface.
    We simply put the command line we would run on a console.
    """
    import nipype.pipeline.engine as pe
    from nipype.interfaces.base import CommandLine

    temporary_env = os.environ.copy()
    temporary_env["SUBJECTS_DIR"] = str(subjects_directory)

    segmentation = pe.Node(
        interface=CommandLine(
            command=f"gtmseg --s {freesurfer_id} --no-seg-stats --xcerseg",
            terminal_output="stream",
        ),
        name="gtmseg",
    )
    segmentation.inputs.environ = temporary_env
    segmentation.run()


def perform_gtmseg(
    caps_dir: Path, subject_id: str, session_id: str, is_longitudinal: bool
):
    """Perform Freesurfer gtmseg.
    It is a command used to perform a segmentation used in some partial volume correction methods.

    Parameters
    ----------
    caps_dir : Path
        CAPS directory
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
    subjects_dir, freesurfer_id = _get_new_subjects_dir(
        is_longitudinal, caps_dir, subject_id, session_id
    )

    freesurfer_mri_folder = subjects_dir / freesurfer_id / "mri"
    gtmseg_file_path = freesurfer_mri_folder / "gtmseg.mgz"

    if not gtmseg_file_path.exists():
        _run_gtmseg(freesurfer_id, subjects_dir)

    # We specify the out file to be in the current directory of execution (easy for us to look at it afterward in the
    # working directory). We copy then the file.
    out_file = Path.cwd() / "gtmseg.mgz"
    copy_file(gtmseg_file_path, out_file)

    # Remove files created during segmentation in the CAPS but not needed
    for filename in ("gtmseg.ctab", "gtmseg.lta"):
        filepath = freesurfer_mri_folder / filename
        if filepath.exists():
            filepath.unlink()
    return out_file


def remove_nan_from_image(image_path: Path) -> Path:
    """Remove NaN values from the provided nifti image.
    This is needed after a registration performed by 'spmregister' : instead
    of filling space with 0, nan are used to extend the PET space.
    We propose to replace them with 0s.
    Parameters
    ----------
    image_path : Path
        The path to the Nifti volume where NaNs need to be replaced by zeros.
    Returns
    -------
    output_image_path : Path
        The path to the volume in Nifti that does not contain any NaNs.
    """
    from clinica.utils.filemanip import get_filename_no_ext

    image = nib.load(image_path)
    data = np.nan_to_num(image.get_fdata(dtype="float32"))
    output_image = nib.Nifti1Image(data, image.affine, header=image.header)
    output_image_path = Path.cwd() / f"no_nan_{get_filename_no_ext(image_path)}.nii.gz"
    nib.save(output_image, output_image_path)
    return output_image_path


def _read_region_source_dst_csv(csv_file: Path) -> pd.DataFrame:
    expected_columns = ["REGION", "SOURCE", "DST"]
    if not csv_file.is_file():
        raise IOError(f"The provided CSV file {csv_file} does not exist.")
    df = pd.read_csv(csv_file, sep=",")
    if df.columns.values.tolist() != expected_columns:
        raise Exception(
            f"CSV file {csv_file} is not in the correct format. "
            f"Columns should be: {expected_columns}."
        )
    return df


def _load_source_dest_region(
    csv_mapping_file: Path,
) -> tuple[np.ndarray, np.ndarray, list[str]]:
    df = _read_region_source_dst_csv(csv_mapping_file)
    return (
        np.asanyarray(list(df.SOURCE)).astype("int"),
        np.asarray(list(df.DST)).astype("int"),
        list(df.REGION),
    )


def _check_mapping_integrity(
    original_labels: np.ndarray, mapping_source_values: np.ndarray
) -> None:
    # Check that each label of original volume (old_label) has a matching transformation in the mapping file
    for label in original_labels:
        if label not in mapping_source_values:
            raise Exception(
                f"Could not find label {label} on conversion table. Add it manually in CSV file to correct error"
            )


def _apply_mapping(
    gtm_segmentation_labels: np.ndarray, source: np.ndarray, dest: np.ndarray
) -> np.ndarray:
    # todo : in test verify the dtype should be integers16
    new_labels_volume = np.zeros(
        gtm_segmentation_labels.shape, dtype=gtm_segmentation_labels.dtype
    )
    # Computing the transformation
    for i, src in enumerate(source):
        new_labels_volume[gtm_segmentation_labels == src] = dest[i]

    return new_labels_volume


def _are_almost_equal(a: float, b: float, rel_tol=1e-9, abs_tol=0.0) -> bool:
    """Measure equality between to floating or double numbers, using 2 thresholds : a
    relative tolerance, and an absolute tolerance
    """
    return abs(a - b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)


def _check_sum(control_image_data: np.ndarray):
    """The sum of a voxel location across the fourth dimension should be 1."""
    sum_voxel_mean = float(sum(sum(sum(control_image_data)))) / control_image_data.size
    if not _are_almost_equal(1.0, sum_voxel_mean):
        raise ValueError(
            "Problem during parcellation: the mean sum of a voxel across "
            f"4th dimension is {sum_voxel_mean} instead of 1.0"
        )


def convert_labels(gtmseg_file: Path, csv_mapping_file: Path) -> list[Path]:
    """Method used on the segmentation from gtmsegmentation.

    The purpose is to reduce the number of label.
    The gathering of labels is specified in a separate file.

    Parameters
    ----------
    gtmseg_file : Path
        The path to the Nifti volume containing the gtmseg segmentation.

    csv_mapping_file : Path
        The path to the mapping .csv file that contains 3 columns : REGION SOURCE DST.
        Separator is , (coma).

    Returns
    -------
    list of Path :
        List of path to the converted volumes according to the .csv file.
        Each volume is a mask representing an area.
    """
    gtm_segmentation = nib.load(gtmseg_file)
    gtm_segmentation.header.set_data_dtype("int8")
    gtm_segmentation_volume = gtm_segmentation.get_fdata(dtype="float32").astype(
        "int16"
    )
    control_volume = np.zeros(gtm_segmentation_volume.shape)

    source, dest, region = _load_source_dest_region(csv_mapping_file)
    _check_mapping_integrity(np.unique(gtm_segmentation_volume), source)

    new_labels_volume = _apply_mapping(gtm_segmentation_volume, source, dest)

    list_of_regions = list()

    # For each label, create a volume file filled with 0s and 1s and save it in current directory under whatever name
    for new_label in np.unique(new_labels_volume):
        region_volume = np.zeros(gtm_segmentation_volume.shape, dtype="uint8")
        region_volume[new_labels_volume == new_label] = 1
        output_path = Path(f"{new_label}.nii.gz").resolve()
        nib.save(
            nib.Nifti1Image(
                region_volume, gtm_segmentation.affine, header=gtm_segmentation.header
            ),
            output_path,
        )
        list_of_regions.append(output_path)
        control_volume = control_volume + region_volume

    _check_sum(control_volume)

    return list_of_regions


def _set_script_for_spm_standalone(
    target_image: Path,
    deformation_field: Path,
    img: Path,
    prefix: str = "subject_space_",
) -> str:
    from textwrap import dedent

    script_file = dedent(
        """
        spm('Defaults', 'fMRI');
        spm_jobman('initcfg');

        jobs{{1}}.spm.util.defs.comp{{1}}.inv.comp{{1}}.def   = {{'{deformation_field}'}};
        jobs{{1}}.spm.util.defs.comp{{1}}.inv.space           = {{'{target}'}};
        jobs{{1}}.spm.util.defs.out{{1}}.pull.fnames          = {{'{img}'}};
        jobs{{1}}.spm.util.defs.out{{1}}.pull.savedir.saveusr = {{'{output_dir}'}};
        jobs{{1}}.spm.util.defs.out{{1}}.pull.interp          = 4;
        jobs{{1}}.spm.util.defs.out{{1}}.pull.mask            = 1;
        jobs{{1}}.spm.util.defs.out{{1}}.pull.fwhm            = [0 0 0];
        jobs{{1}}.spm.util.defs.out{{1}}.pull.prefix          = '{prefix}';

        spm_jobman('run', jobs);
        """
    )

    script_file = script_file.format(
        deformation_field=deformation_field,
        target=target_image,
        img=img,
        output_dir=str(Path.cwd()),
        prefix=prefix,
    )

    return script_file


def _call_spm_standalone(script_location: Path, expected_output_location: Path) -> str:
    from clinica.utils.check_dependency import get_spm_standalone_home
    from clinica.utils.spm import _get_real_spm_standalone_file

    # TODO : This might not even be needed with cmd line setting done before in the pipeline
    spm_file = _get_real_spm_standalone_file(get_spm_standalone_home())
    cmdline = f"$SPMSTANDALONE_HOME/{spm_file} $MCR_HOME batch {script_location}"

    run_command_as_subprocess(
        "runApplyInverseDeformationField_SPM_standalone",
        cmdline,
        expected_output_location,
    )

    return cmdline


def run_apply_inverse_deformation_field_SPM_standalone(
    target_image: Path, deformation_field: Path, img: Path
) -> Path:
    """
    We directly create a batch file that SPM standalone can run. This function does not check whether SPM standalone must be used. Previous
    check when building the pipeline ensures that all the env vars exists ($SPMSTANDALONE_HOME and $MCR_HOME)

    Parameters
    ----------
    target_image : Path
        Path to the target image
    deformation_field : Path
        Path to the deformation field
    img : Path
        Path to the moving image, ie which is warped

    Returns
    -------
    Path :
        Path to the result
    """
    prefix = "subject_space_"

    # Write SPM batch command directly in a script that is readable by SPM standalone
    script_location = Path("./m_script.m").resolve()
    script_file = _set_script_for_spm_standalone(
        target_image, deformation_field, img, prefix
    )

    with open(script_location, "w", encoding="utf-8") as f:
        f.write(script_file)

    output_file = (
        Path.cwd() / f"{prefix}{img.name}"
    )  # TODO : if issue with symlinks use .resolve()

    _call_spm_standalone(script_location, output_file)

    return output_file


def _get_eroded_mask_and_size(mask: Path) -> tuple[np.ndarray, int]:
    from clinica.utils.exceptions import ClinicaImageError

    eroded_mask_nifti = nib.load(mask)
    eroded_mask = eroded_mask_nifti.get_fdata(dtype="float32")
    eroded_mask = eroded_mask > 0

    if (mask_size := np.sum(eroded_mask)) == 0:
        raise ClinicaImageError(
            f"The eroded mask located at {mask} contains only zero values. "
            "A problem likely occurred when moving the eroded mask from MNI to gtmsegspace."
        )
    return eroded_mask, mask_size


def _get_mean_pet_activity_within_mask(mask: Path, pet_data: np.ndarray) -> float:
    eroded_mask, mask_size = _get_eroded_mask_and_size(mask)
    return np.sum(eroded_mask * pet_data) / mask_size


def normalize_suvr(pet_path: Path, mask: Path) -> Path:
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

    Returns
    -------
    Path :
        The path to the SUVR normalized volume in the current directory.

    Raises
    ------
    ClinicaImageError :
        If the provided eroded mask contains only zero values.
    """
    # Load PET data (they must be in gtmsegspace, or same space as label file)
    pet_image_nifti = nib.load(pet_path)
    pet_data = pet_image_nifti.get_fdata(dtype="float32")

    # Mask unwanted values to determine mean uptake value
    mean_pons_pet_activity = _get_mean_pet_activity_within_mask(mask, pet_data)

    # Then normalize PET data by this mean activity
    suvr_image_nifti = nib.Nifti1Image(
        pet_data / mean_pons_pet_activity,
        pet_image_nifti.affine,
        header=pet_image_nifti.header,
    )
    suvr_filename = Path.cwd() / f"suvr_{pet_path.name}"
    nib.save(suvr_image_nifti, suvr_filename)
    return suvr_filename


def _make_freesurfer_command_mac_compatible(command: str) -> str:
    return "export DYLD_LIBRARY_PATH=$FREESURFER_HOME/lib/gcc/lib && " + command


def _build_mris_expand_cmd(in_surface: Path) -> str:
    cmd = f"mris_expand -thickness -N 13 {in_surface} 0.65 {in_surface.name}_exp-"
    # If system is MacOS, this export command must be run just before the mri_vol2surf command to bypass MacOs security
    if platform.system().lower().startswith("darwin"):
        cmd = _make_freesurfer_command_mac_compatible(cmd)

    return cmd


def _get_numbered_exp_filename(filename: Path, number: int) -> Path:
    # Expects a filename like Path.cwd() / lh.white to output Path.cwd() / lh.white_exp-00N
    return filename.with_name(f"{filename.name}_exp-{str(number).zfill(3)}")


def _check_mri_expand_file_location_then_move(
    working_directory: Path, input_file_location: Path
) -> Path:
    expected_location = working_directory / input_file_location.name
    if _get_numbered_exp_filename(input_file_location, 0).is_file():
        for n in range(0, 14):
            _get_numbered_exp_filename(input_file_location, n)
            move_file(
                _get_numbered_exp_filename(input_file_location, n),
                _get_numbered_exp_filename(expected_location, n),
            )

    return expected_location


def run_mris_expand(surface: Path) -> list[Path]:
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

    Returns
    -------
    List of Path :
        List of path to the generated surfaces.

    Notes
    -----
    'mris_expand' write results where the script is executed

    -N is a hidden parameter (not documented) that allows the user to specify
    the number of surface generated between source and final target surface.
    Here target is 65% of thickness, with 13 surfaces.
    Then we only keep the surfaces we are interested in.
    """

    run_command_as_subprocess("mris_expand", _build_mris_expand_cmd(surface))

    # Remove useless surfaces (0%, 5%, 10%, 15%, 20%, 25% and 30% of thickness)
    cprint(msg="Removing unnecessary mris_expands outputs (000 to 007)", lvl="debug")

    expected_location = _check_mri_expand_file_location_then_move(
        working_directory=Path.cwd(), input_file_location=surface
    )

    for file in [_get_numbered_exp_filename(expected_location, x) for x in range(0, 7)]:
        file.unlink()

    return [_get_numbered_exp_filename(expected_location, x) for x in range(7, 14)]


def _build_mri_surf2surf_command(
    surface: Path,
    registration: Path,
    gtmsegfile: Path,
    freesurfer_id: str,
    output_file: Path,
) -> str:
    hemisphere = HemiSphere(surface.name[0:2])
    surface_name = surface.name[3:]
    command = f"mri_surf2surf --reg {registration} {gtmsegfile} --sval-xyz {surface_name} --hemi {hemisphere.value} --tval-xyz {gtmsegfile} --tval {output_file} --s {freesurfer_id}"

    # If system is MacOS, this export command must be run just
    # before the mri_vol2surf command to bypass MacOs security.
    if platform.system().lower().startswith("darwin"):
        command = _make_freesurfer_command_mac_compatible(command)
    return command


def run_mri_surf2surf(
    surface: Path,
    registration: Path,
    gtmsegfile: Path,
    subject_id: str,
    session_id: str,
    caps_dir: Path,
    is_longitudinal: bool,
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

    Returns
    -------
    Path :
        The path to the converted surface in current directory.
    """

    subject_directory, freesurfer_id = _get_new_subjects_dir(
        is_longitudinal, caps_dir, subject_id, session_id
    )
    copy_file(surface, subject_directory / freesurfer_id / "surf", exist_ok=True)
    output_path = Path.cwd() / f"{surface.name}_gtmsegspace"

    temporary_env = os.environ.copy()
    temporary_env["SUBJECTS_DIR"] = str(subject_directory)

    run_command_as_subprocess(
        "mri_surf2surf",
        _build_mri_surf2surf_command(
            surface, registration, gtmsegfile, freesurfer_id, output_path
        ),
        env=temporary_env,
    )

    (subject_directory / freesurfer_id / "surf" / surface.name).unlink(missing_ok=False)

    return output_path


def _build_mri_vol2surf_command(
    pet_volume: Path,
    surface: Path,
    freesurfer_id: str,
    output_file: Path,
) -> str:
    hemisphere = HemiSphere(surface.name[0:2])
    surface_name = surface.name[3:]
    command = (
        f"mri_vol2surf --mov {pet_volume} --o {output_file} --surf {surface_name} --hemi {hemisphere.value} "
        f"--regheader {freesurfer_id} --ref gtmseg.mgz --interp nearest"
    )

    if platform.system().lower().startswith("darwin"):
        command = _make_freesurfer_command_mac_compatible(command)

    return command


def run_mri_vol2surf(
    pet_volume: Path,
    surface: Path,
    subject_id: str,
    session_id: str,
    caps_dir: Path,
    gtmsegfile: Path,
    is_longitudinal: bool,
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

    Returns
    -------
    Path :
        The path to the data projected onto the surface.
    """

    subjects_dir, freesurfer_id = _get_new_subjects_dir(
        is_longitudinal, caps_dir, subject_id, session_id
    )

    copy_file(surface, subjects_dir / freesurfer_id / "surf")
    gtmsegfile_copy = subjects_dir / freesurfer_id / "mri" / "gtmseg.mgz"
    if not gtmsegfile_copy.exists():
        copy_file(gtmsegfile, gtmsegfile_copy)

    # TODO write nicer way to grab hemi (difficulty caused by the dots in filenames)
    # extract hemisphere based on filename
    hemisphere = HemiSphere(surface.name[0:2])
    output_file = Path.cwd() / f"{hemisphere.value}.projection_{surface.name}.mgh"

    temporary_env = os.environ.copy()
    temporary_env["SUBJECTS_DIR"] = str(subjects_dir)

    run_command_as_subprocess(
        "mri_vol2surf",
        _build_mri_vol2surf_command(pet_volume, surface, freesurfer_id, output_file),
        env=temporary_env,
    )

    (subjects_dir / freesurfer_id / "surf" / surface.name).unlink(missing_ok=False)
    # TODO careful here...
    # Removing gtmseg.mgz may lead to problems as other vol2surf are using it
    gtmsegfile_copy.unlink(missing_ok=False)

    return output_file


def _build_weighted_mean_surface_image(surfaces: Sequence[Path]) -> nib.MGHImage:
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
    tuple[float, float, float, float, float, float, float]
):
    """TODO: Find out where do these numbers come from.."""
    return 0.1034, 0.1399, 0.1677, 0.1782, 0.1677, 0.1399, 0.1034


def compute_weighted_mean_surface(surfaces: Sequence[Path]) -> Path:
    """Compute a weighted average at each node of the surface.

    The weight are defined by a normal distribution (centered on the mid surface).

    Parameters
    ----------
    surfaces : Sequence of Path
        The paths to the data projected on the 7 surfaces (35 to 65 % of thickness) at each nodes.

    Returns
    -------
    Path :
        The path to the data averaged.
    """
    _assert_seven_surfaces(surfaces)
    out_surface = (
        Path.cwd()
        / f"{HemiSphere(surfaces[0].name[0:2]).value}.averaged_projection_on_cortical_surface.mgh"
    )
    nib.save(_build_weighted_mean_surface_image(surfaces), out_surface)

    return out_surface


def _fsaverage_was_copied(src: Path, dst: Path) -> bool:
    # todo : does it ever happen that it was not copied bc already exists ?
    if not (dst / "fsaverage").exists():
        shutil.copytree(src / "fsaverage", dst / "fsaverage")
        return True
    return False


def project_onto_fsaverage(
    projection: Path,
    subject_id: str,
    session_id: str,
    caps_dir: Path,
    fwhm: int,
    is_longitudinal: bool,
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

    Returns
    -------
    Path :
        The path to the data averaged.
    """
    subjects_dir_backup = os.path.expandvars("$SUBJECTS_DIR")

    subjects_dir, freesurfer_id = _get_new_subjects_dir(
        is_longitudinal, caps_dir, subject_id, session_id
    )

    # copy fsaverage folder next to : subject_id + '_' + session_id
    # for the mris_preproc command to properly find src and target
    fsaverage_has_been_copied = _fsaverage_was_copied(
        Path(subjects_dir_backup), subjects_dir
    )

    # also copy the mgh file in the surf folder (needed by MRISPreproc)
    projection_in_surf_folder = subjects_dir / freesurfer_id / "surf" / projection.name

    if not projection_in_surf_folder.exists():
        copy_file(projection, projection_in_surf_folder)

    out_fsaverage = Path.cwd() / f"fsaverage_fwhm-{fwhm}_{projection.name}"

    _run_mris_preproc_as_standalone_nipype_node(
        projection, freesurfer_id, fwhm, out_fsaverage, subjects_directory=subjects_dir
    )

    # remove projection file from surf folder
    projection_in_surf_folder.unlink(missing_ok=False)

    # remove fsaverage if it has been copied
    if fsaverage_has_been_copied:
        shutil.rmtree(subjects_dir / "fsaverage")

    return out_fsaverage


def _run_mris_preproc_as_standalone_nipype_node(
    projection: Path,
    freesurfer_id: str,
    fwhm: float,
    output_file: Path,
    subjects_directory: Path,
):
    from nipype.interfaces.freesurfer import MRISPreproc

    temporary_env = os.environ.copy()
    temporary_env["SUBJECTS_DIR"] = str(subjects_directory)

    projection_node = MRISPreproc()
    projection_node.inputs.target = "fsaverage"
    projection_node.inputs.subjects = [freesurfer_id]
    projection_node.inputs.fwhm = fwhm
    projection_node.inputs.hemi = HemiSphere(projection.name[0:2]).value
    projection_node.inputs.surf_measure = projection.name[3:]
    projection_node.inputs.out_file = str(output_file)
    projection_node._environ = temporary_env
    projection_node.run()


def _assert_seven_surfaces(surfaces: Sequence[Path]):
    if (n_surfaces := len(surfaces)) != 7:
        raise ValueError(
            "There should be 7 surfaces at this point of the pipeline. "
            f"However 'compute_weighted_mean_surface' received {n_surfaces} surfaces. "
            "Something probably went wrong in prior steps of the pipeline."
        )


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


def reformat_surfname(
    hemisphere: HemiSphere, left_surface: Path, right_surface: Path
) -> Path:
    if hemisphere == HemiSphere.LEFT:
        return left_surface
    if hemisphere == HemiSphere.RIGHT:
        return right_surface


def compute_average_pet_signal_based_on_annotations(
    pet_projections: tuple[Path, Path], atlas_files: tuple[str, FreeSurferAnnotation]
) -> list[Path]:
    """Computes the average of PET signal based on annot files from Freesurfer.

    Those files describe the brain according to known atlases.

    Parameters
    ----------
    pet_projections : tuple of two Path
        The paths to the PET projection (must be a MGH file) [left_hemisphere, right_hemisphere].

    atlas_files : tuple[str, FreeSurferAnnotationImage]
        Tuple containing path to lh and rh annotation files for any number of atlases and their names.

    Returns
    -------
    Path :
        The path to the tsv containing average PET values.

    Raises
    ------
    ValueError :
        If not exactly two files were provided through the argument 'pet_projections'.
    """
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

    for annotation in atlas_files:
        annotation.replace_minus_one_annotation_with_zero()
        average_region = []
        for region_id, reg_name in enumerate(annotation.region_names):
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
                "mean_scalar": list(average_region),
            }
        )

        filename_atlas_tsv = Path.cwd() / f"{annotation.atlas_name}.tsv"
        filename_tsv.append(filename_atlas_tsv)
        final_tsv.to_csv(
            filename_atlas_tsv,
            sep="\t",
            index=False,
            columns=["index", "label_name", "mean_scalar"],
        )
    return filename_tsv


def merge_nifti_volumes(inputs: list[str]) -> str:
    import os

    import nibabel as nib
    from nilearn.image import concat_imgs

    sorted_inputs = sorted(
        inputs, key=lambda p: int(p.split("/")[-1].split(".nii.gz")[0])
    )
    merged_image = concat_imgs([nib.load(p) for p in sorted_inputs])
    output_path = os.getcwd() + "/merged_image.nii.gz"
    nib.save(merged_image, output_path)
    return output_path


def get_regexp_substitutions(
    pet_tracer: Tracer,
    region: SUVRReferenceRegion,
    is_longitudinal: bool,
) -> list[tuple[str, str]]:
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


def _get_mid_surface_substitutions(is_longitudinal: bool) -> tuple[str, str]:
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
) -> tuple[str, str]:
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
) -> tuple[str, str]:
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
) -> tuple[str, str]:
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
