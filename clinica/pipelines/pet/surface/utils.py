from os import PathLike
from pathlib import Path
from typing import Union

__all__ = [
    "get_new_subjects_dir",
    "perform_gtmseg",
    "remove_nan",
    "make_label_conversion",
    "run_mri_vol2surf",
    "compute_weighted_mean_surface",
    "project_onto_fsaverage",
    "get_mid_surface",
    "reformat_surfname",
    "produce_tsv",
    "merge_nifti_volumes",
    "run_ApplyInverseDeformationField_SPM_standalone",
    "run_mri_surf2surf",
    "normalize_suvr",
    "mris_expand",
    "remove_nan_from_image",
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


def get_new_subjects_dir(
    is_longitudinal: bool,
    caps_dir: Union[str, PathLike],
    subject_id: str,
    session_id: str,
):
    """Extract SUBJECT_DIR.

    Extract path to FreeSurfer segmentation in CAPS folder and FreeSurfer ID
    (e.g. sub-CLNC01_ses-M000.long.sub-CLNC01_long-M000M018 or sub-CLNC01_ses-M000).
    """
    caps_dir = Path(caps_dir)
    # todo : unclear what this function is supposed to do
    # todo : t1 or pet ?

    root = caps_dir / "subjects" / subject_id / session_id / "t1"

    if is_longitudinal:
        long_folds = _get_longitudinal_folder_name(root)

        return (
            root / long_folds / "freesurfer_longitudinal",
            f"{subject_id}_{session_id}.long.{subject_id}_{long_folds}",
        )
    return root / "freesurfer_cross_sectional", subject_id + "_" + session_id


def perform_gtmseg(caps_dir, subject_id, session_id, is_longitudinal):
    """gtmseg is a freesurfer command used to perform a segmentation used in some partial volume correction methods.

    Warnings:
        - This method changes the environment variable $SUBJECTS_DIR (but put
          the original one back after execution).  This has not been intensely
          tested whether it can lead to some problems : (for instance if 2
          subjects are running in parallel)

    Args:
        (string) caps_dir : CAPS directory.
        (string) subject_id: The subject_id (something like sub-ADNI002S4213)
        (string) session_id: The session id ( something like : ses-M012)
        (bool)   is_longitudinal: If longitudinal processing, subjects_dir must be put elsewhere

    Returns:
        (string) Path to the segmentation volume : a volume where each voxel
        has a label (ranging [0 2035] ), see Freesurfer lookup table to see the
        labels with their corresponding names.
    """
    import os
    import shutil

    import nipype.pipeline.engine as pe
    from nipype.interfaces.base import CommandLine

    # Old subject_dir is saved for later
    subjects_dir_backup = os.path.expandvars("$SUBJECTS_DIR")

    root_env, freesurfer_id = get_new_subjects_dir(
        is_longitudinal, caps_dir, subject_id, session_id
    )

    # Set the new subject dir for the function to work properly
    os.environ["SUBJECTS_DIR"] = root_env

    if not os.path.exists(
        os.path.join(
            os.path.expandvars("$SUBJECTS_DIR"), freesurfer_id, "mri", "gtmseg.mgz"
        )
    ):
        # Creation of standalone node based on Command Line Interface.
        # We simply put the command line we would run on a console
        segmentation = pe.Node(
            interface=CommandLine(
                "gtmseg --s " + freesurfer_id + " --no-seg-stats --xcerseg",
                terminal_output="stream",
            ),
            name="gtmseg",
        )
        segmentation.run()

    # We specify the out file to be in the current directory of execution (easy for us to look at it afterward in the
    # working directory). We copy then the file.
    out_file = os.path.abspath("./gtmseg.mgz")
    shutil.copy(
        os.path.join(
            os.path.expandvars("$SUBJECTS_DIR"), freesurfer_id, "mri", "gtmseg.mgz"
        ),
        out_file,
    )

    # Remove bunch of files created during segmentation in caps dir and not needed
    gtmsegcab = os.path.join(
        os.path.expandvars("$SUBJECTS_DIR"), freesurfer_id, "mri", "gtmseg.ctab"
    )
    if os.path.exists(gtmsegcab):
        os.remove(gtmsegcab)

    gtmseglta = os.path.join(
        os.path.expandvars("$SUBJECTS_DIR"), freesurfer_id, "mri", "gtmseg.lta"
    )
    if os.path.exists(gtmseglta):
        os.remove(gtmseglta)

    # Set back the SUBJECT_DIR environment variable of the user
    os.environ["SUBJECTS_DIR"] = subjects_dir_backup
    return out_file


def remove_nan(volname):
    """remove_nan is a method needed after a registration performed by spmregister : instead of filling space with 0, nan
    are used to extend the PET space. We propose to replace them with 0s.

    Args:
        (string) volname : path to the Nifti volume where NaNs need to be replaced by 0s

    Returns:
        (string) Path to the volume in Nifti that does not contain any NaNs
    """

    # todo : discriminate with remove_nan_from_image
    import os

    import nibabel as nib
    import numpy as np

    # Load the volume and get the data
    nifti_in = nib.load(volname)
    data = np.nan_to_num(nifti_in.get_fdata(dtype="float32"))

    # Now create final image (using header of original image), and save it in current directory
    nifti_out = nib.Nifti1Image(data, nifti_in.affine, header=nifti_in.header)
    filename = os.path.basename(volname)
    vol_wo_nan = "./no_nan_" + filename + ".gz"
    vol_wo_nan = os.path.abspath(vol_wo_nan)
    nib.save(nifti_out, vol_wo_nan)
    return vol_wo_nan


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


def make_label_conversion(gtmsegfile, csv):
    """make_label_conversion is a method used on the segmentation from gtmsegmentation. The purpose is to reduce the
    number of label. The gathering of labels is specified in a separate file

    Args:
        (string) gtmsegfile   : path to the Nifti volume containing the gtmseg segmentation
        (string) csv          : path to .csv file that contains 3 columns : REGION SOURCE DST. Separator is , (coma).

    Returns:
        (list of strings) List of path to the converted volumes according to the .csv file. Each volume is a mask
        representing an area
    """
    import os

    import nibabel as nib
    import numpy
    import pandas

    def isclose(a, b, rel_tol=1e-9, abs_tol=0.0):
        """Small function designed to measure equality between to floating or double numbers, using 2 thresholds : a
        relative tolerance, and an absolute tolerance
        """
        return abs(a - b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

    # Read label from gtmsegfile, change data into integers in order to have no problems when testing equality of labels
    label = nib.load(gtmsegfile)
    label.header.set_data_dtype("int8")
    volume = label.get_fdata(dtype="float32")

    # Unique function gives a list where each label of the volume is listed once
    old_labels = numpy.unique(volume)
    old_labels = old_labels.astype("int16")

    # allsum is a control volume (sum of a pixel across the 4 th dimension must be equal to 1)
    allsum = numpy.zeros(volume.shape)

    # Reading of csv file, raise exception if the pattern REGION, SOURCE, DST is not found
    if not os.path.isfile(csv):
        raise Exception("The CSV file does not exist.")
    convert_lut = pandas.io.parsers.read_csv(csv, sep=",")
    if list(convert_lut.columns.values) != ["REGION", "SOURCE", "DST"]:
        raise Exception(
            f"CSV file {csv} is not in the correct format. Columns should be: REGION, SOURCE, DST"
        )

    # Extract columns to a list form (values converted into integers)
    src = list(convert_lut.SOURCE)
    src_val = numpy.asanyarray(src)
    src_val = src_val.astype("int")

    dst = list(convert_lut.DST)
    dst_val = numpy.asarray(dst)
    dst_val = dst_val.astype("int")

    # Check that each label of original volume (old_label) has a matching transformation in the csv file
    for i in range(old_labels.size):
        index = numpy.argwhere(src_val == old_labels[i])
        # Size 0 means no occurrence found
        if index.size == 0:
            raise Exception(
                f"Could not find label {old_labels[i]} on conversion table. Add it manually in CSV file to correct error"
            )

    # Instantiation of final volume, with same dtype as original volume
    new_volume = numpy.zeros(volume.shape, dtype=volume.dtype)
    # Computing the transformation
    for i in range(len(src)):
        new_volume[volume == src_val[i]] = dst_val[i]
    # Get unique list of new label
    new_labels = numpy.unique(new_volume)
    new_labels = new_labels.astype("int")
    list_of_regions = list()

    # For each label, create a volume file filled with 0s and 1s and save it in current directory under whatever name
    for i in range(new_labels.size):
        region_volume = numpy.zeros(volume.shape, dtype="uint8")
        region_volume[new_volume == new_labels[i]] = 1
        myNifti = nib.Nifti1Image(region_volume, label.affine, header=label.header)
        current_path = "./" + str(new_labels[i]) + ".nii.gz"
        current_path = os.path.abspath(current_path)
        nib.save(myNifti, current_path)
        list_of_regions.append(current_path)
        allsum = allsum + region_volume

    # The sum of a voxel location across the fourth dimension should be 1
    sum_voxel_mean = float(sum(sum(sum(allsum)))) / allsum.size
    if not isclose(1.0, sum_voxel_mean):
        raise Exception(
            f"Problem during parcellation: the mean sum of a voxel across 4th dimension is {sum_voxel_mean} instead of 1.0"
        )
    # The list of files is returned
    return list_of_regions


def run_ApplyInverseDeformationField_SPM_standalone(target, deformation_field, img):
    """
    We directly create a batch file that SPM standalone can run. This function does not check whether SPM standalone must be used. Previous
    check when building the pipeline ensures that all the env vars exists ($SPMSTANDALONE_HOME and $MCR_HOME)
    """
    import os
    import subprocess
    from os.path import abspath, basename, exists, join
    from textwrap import dedent

    from clinica.utils.check_dependency import get_spm_standalone_home
    from clinica.utils.spm import _get_real_spm_standalone_file

    prefix = "subject_space_"

    # Write SPM batch command directly in a script that is readable by SPM standalone
    script_location = abspath("./m_script.m")
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
        target=target,
        img=img,
        output_dir=abspath(os.getcwd()),
        prefix=prefix,
    )

    with open(script_location, "w", encoding="utf-8") as f:
        f.write(script_file)

    # TODO : This might not even be needed with cmd line setting done prior
    spm_file = _get_real_spm_standalone_file(get_spm_standalone_home())
    cmdline = f"$SPMSTANDALONE_HOME/{spm_file} $MCR_HOME batch {script_location}"

    subprocess_run = subprocess.run(
        cmdline,
        shell=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    if code := subprocess_run.returncode != 0:
        raise ValueError(
            f"runApplyInverseDeformationField_SPM_standalone failed, returned non-zero code with {code}"
        )

    output_file = join(abspath("./"), prefix + basename(img))
    if not exists(output_file):
        raise IOError(
            "Something went wrong while trying to run runApplyInverseDeformationField_SPM_standalone"
            + ". Output file not generated. Command launched :\n\t "
            + cmdline
            + "\n. We strongly recommend that you use the supported version of Matlab MCR "
            + " recommended by the creators of SPM."
        )
    return output_file


def normalize_suvr(pet_path, mask):
    """normalize_suvr is a way of getting suvr from your pet image, based on the segmentation performed by
    gtmsegmentation. The Standard Uptake Value ratio is computed by dividing the whole PET volume by the mean value
    observed in the pons.

    Args:
        (string) pet_path     : path to the Nifti volume containing PET scan, realigned on upsampled T1
        (string) mask         : mask of the pons (18FFDG) or pons+cerebellum (18FAV45) already eroded

    Returns:
        (string) Path to the suvr normalized volume in the current directory
    """
    import os

    import nibabel as nib

    # Load mask
    eroded_mask_nifti = nib.load(mask)
    eroded_mask = eroded_mask_nifti.get_fdata(dtype="float32")
    eroded_mask = eroded_mask > 0

    # Load PET data (they must be in gtmsegspace, or same space as label file)
    pet = nib.load(pet_path)
    pet_data = pet.get_fdata(dtype="float32")

    # check that eroded mask is not null
    mask_size = sum(sum(sum(eroded_mask)))
    if mask_size == 0:
        raise Exception(
            "Number of non-zero value of mask is 0. A problem occurred when moving the eroded mask from MNI to gtmsegspace"
        )

    # Mask unwanted values to determine mean uptake value
    pons_pet_activity = eroded_mask * pet_data
    mean_pons_pet_activity = sum(sum(sum(pons_pet_activity))) / mask_size

    # Then normalize PET data by this mean activity
    suvr_pet_data = pet_data / mean_pons_pet_activity
    suvr = nib.Nifti1Image(suvr_pet_data, pet.affine, header=pet.header)
    suvr_filename = "suvr_" + os.path.basename(pet_path)
    suvr_filename = os.path.abspath("./" + suvr_filename)
    nib.save(suvr, suvr_filename)
    return suvr_filename


def _setting_mris_expand_cmd(in_surface) -> str:
    from pathlib import Path
    from sys import platform

    cmd = (
        "mris_expand -thickness -N 13 "
        + in_surface
        + " 0.65 "
        + Path(in_surface).name
        + "_exp-"
    )
    # If system is MacOS, this export command must be run just before the mri_vol2surf command to bypass MacOs security
    if platform == "darwin":
        cmd = "export DYLD_LIBRARY_PATH=$FREESURFER_HOME/lib/gcc/lib && " + cmd

    return cmd


def _running_mris_expand_with_subprocess(cmd: str) -> None:
    import subprocess

    subprocess_mris_expand = subprocess.run(
        cmd,
        shell=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    if subprocess_mris_expand.returncode != 0:
        raise ValueError("mris_expand failed, returned non-zero code")


def _check_mri_expand_file_location_then_move(
    working_directory: str, input_file_location: str
) -> str:
    import shutil
    from pathlib import Path

    filename = Path(input_file_location).name
    expected_location = f"{working_directory}/{filename}_exp-"

    if Path(input_file_location + "_exp-000").is_file():
        for i in range(0, 14):
            identifier = str(i).zfill(3)
            shutil.move(
                f"{input_file_location}_exp-" + identifier,
                expected_location + identifier,
            )

    return expected_location


def mris_expand(in_surface):
    """mris_expand is using the freesurfer function of the same name. It expands the white input surface toward the pial,
    generating 7 surfaces at 35%, 40%, 45%, 50%, 55%, 60%, 65% of thickness.

    Args:
        (string) in_surface : Path to the input white surface, but must be named lh.white or rh.white, and the folder
            containing the surface file must also have ?h.pial, ?.sphere, ?h.thickness (freesurfer surf folder)

    Returns:
        (list of strings) List of path to the generated surfaces
    """
    import os

    from pipelines.pet.surface.utils import (  # noqa
        _check_mri_expand_file_location_then_move,
        _running_mris_expand_with_subprocess,
        _setting_mris_expand_cmd,
    )

    from clinica.utils.stream import cprint

    # RQ 1 : mris_expand write results where the script is executed
    # RQ 2 : -N is a hidden parameter (not documented) that allows the user to specify the number of surface generated between
    # source and final target surface. Here target is 65% of thickness, with 13 surfaces. Then we only keep the surfaces
    # we are interested in.

    _running_mris_expand_with_subprocess(_setting_mris_expand_cmd(in_surface))

    # Remove useless surfaces (0%, 5%, 10%, 15%, 20%, 25% and 30% of thickness)
    cprint(msg="Removing unnecessary mris_expands outputs (000 to 007)", lvl="debug")

    out_file = _check_mri_expand_file_location_then_move(
        working_directory=os.getcwd(), input_file_location=in_surface
    )

    for file in [out_file + str(x).zfill(3) for x in range(0, 7)]:
        os.remove(file)

    return [os.path.abspath(out_file + str(x).zfill(3)) for x in range(7, 14)]


def run_mri_surf2surf(
    in_surface, reg_file, gtmsegfile, subject_id, session_id, caps_dir, is_longitudinal
):
    """surf2surf is a wrapper of freesurfer command mri_surf2surf. Here the aim is to convert a input surface (which is
    the native space of the subject), into the same surface but in the gtmseg space (space of the volume generated by
    the gtmsegmentation)

    Args:
        (string) in_surface : surface file that needs to be converted
        (string) reg_file   : Path to a registration file that represents the transformation needed to go from the native
            space to the gtmsegspace (see https://surfer.nmr.mgh.harvard.edu/fswiki/FsAnat-to-NativeAnat for more
            details)
        (string) gtmsegfile : Path to the gtm segmentation file
        (string) subject_id : The subject_id (something like sub-ADNI002S4213)
        (string) session_id : The session id ( something like : ses-M012)
        (string) caps_dir   : Path to the CAPS directory
        (bool)   is_longitudinal: longitudinal files

    Returns:
        (string) Path to the converted surface in current directory
    """
    import os
    import shutil
    import subprocess
    import sys

    from pipelines.pet.surface.utils import get_new_subjects_dir

    # set subjects_dir env. variable for mri_surf2surf to work properly
    subjects_dir_backup = os.path.expandvars("$SUBJECTS_DIR")

    root_env, freesurfer_id = get_new_subjects_dir(
        is_longitudinal, caps_dir, subject_id, session_id
    )

    os.environ["SUBJECTS_DIR"] = root_env

    # make a copy of surface file to surface directory in CAPS in order to allow processing
    shutil.copy(
        in_surface,
        os.path.join(os.path.expandvars("$SUBJECTS_DIR"), freesurfer_id, "surf"),
    )

    # TODO write nicer way to grab hemi & filename (difficulty caused by the dots in filenames)
    # extract hemisphere based on filename
    hemi = os.path.basename(in_surface)[0:2]
    surfname = os.path.basename(in_surface)[3:]

    # Perform surf2surf algorithm
    tval = os.path.abspath("./" + os.path.basename(in_surface) + "_gtmsegspace")
    cmd = (
        "mri_surf2surf --reg %s %s --sval-xyz %s --hemi %s --tval-xyz %s --tval %s --s %s "
        % (reg_file, gtmsegfile, surfname, hemi, gtmsegfile, tval, freesurfer_id)
    )

    # If system is MacOS, this export command must be run just before the mri_vol2surf command to bypass MacOs security
    if sys.platform == "darwin":
        cmd = "export DYLD_LIBRARY_PATH=$FREESURFER_HOME/lib/gcc/lib && " + cmd
    subprocess_mri_surf2surf = subprocess.run(
        cmd,
        shell=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    if subprocess_mri_surf2surf.returncode != 0:
        raise ValueError("mri_surf2surf failed, returned non-zero code")

    # remove file in caps
    os.remove(
        os.path.join(
            os.path.expandvars("$SUBJECTS_DIR"),
            freesurfer_id,
            "surf",
            os.path.basename(in_surface),
        )
    )

    # put back original subjects_dir env
    os.environ["SUBJECTS_DIR"] = subjects_dir_backup

    return tval


def run_mri_vol2surf(
    volume, surface, subject_id, session_id, caps_dir, gtmsegfile, is_longitudinal
):
    """vol2surf is a wrapper of freesurfer command mri_vol2surf. It projects the volume into the surface : the value at
    each vertex is given by the value of the voxel it intersects

    Args:
        (string) volume     : Path to PET volume (in gtmseg space) that needs to be mapped into surface
        (string) surface    : Path to surface file
        (string) gtmsegfile :l Path to the gtm segmentation file (provides information on space, labels are not used
        (string) subject_id : The subject_id (something like sub-ADNI002S4213)
        (string) session_id : The session id ( something like : ses-M012)
        (string) caps_dir   : Path to the CAPS directory

    Returns:
        (string) Path to the data projected onto the surface
    """
    import os
    import shutil
    import subprocess
    import sys

    from pipelines.pet.surface.utils import get_new_subjects_dir

    # set subjects_dir env. variable for mri_vol2surf to work properly
    subjects_dir_backup = os.path.expandvars("$SUBJECTS_DIR")

    root_env, freesurfer_id = get_new_subjects_dir(
        is_longitudinal, caps_dir, subject_id, session_id
    )

    os.environ["SUBJECTS_DIR"] = root_env

    # TODO write nicer way to grab hemi & filename (difficulty caused by the dots in filenames)
    # extract hemisphere based on filename
    hemi = os.path.basename(surface)[0:2]
    surfname = os.path.basename(surface)[3:]

    # copy surface file in caps surf folder to allow processing
    shutil.copy(
        surface,
        os.path.join(os.path.expandvars("$SUBJECTS_DIR"), freesurfer_id, "surf"),
    )

    if not os.path.exists(
        os.path.join(
            os.path.expandvars("$SUBJECTS_DIR"), freesurfer_id, "mri", "gtmseg.mgz"
        )
    ):
        shutil.copy(
            gtmsegfile,
            os.path.join(
                os.path.expandvars("$SUBJECTS_DIR"), freesurfer_id, "mri", "gtmseg.mgz"
            ),
        )

    # execute vol2surf
    output = os.path.abspath(
        "./" + hemi + ".projection_" + os.path.basename(surface) + ".mgh"
    )
    cmd = "mri_vol2surf"
    cmd += " --mov " + volume
    cmd += " --o " + output
    cmd += " --surf " + surfname
    cmd += " --hemi " + hemi
    cmd += " --regheader " + freesurfer_id
    cmd += " --ref gtmseg.mgz"
    cmd += " --interp nearest"

    # If system is MacOS, this export command must be run just before the mri_vol2surf command to bypass MacOs security
    if sys.platform == "darwin":
        cmd = "export DYLD_LIBRARY_PATH=$FREESURFER_HOME/lib/gcc/lib && " + cmd
    subprocess_mri_vol2surf = subprocess.run(
        cmd,
        shell=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    if subprocess_mri_vol2surf.returncode != 0:
        raise ValueError("mri_vol2surf failed, returned non-zero code")

    # remove file in caps
    os.remove(
        os.path.join(
            os.path.expandvars("$SUBJECTS_DIR"),
            freesurfer_id,
            "surf",
            os.path.basename(surface),
        )
    )
    # TODO careful here...
    # Removing gtmseg.mgz may lead to problems as other vol2surf are using it
    os.remove(
        os.path.join(
            os.path.expandvars("$SUBJECTS_DIR"), freesurfer_id, "mri", "gtmseg.mgz"
        )
    )

    # put back original subjects_dir env
    os.environ["SUBJECTS_DIR"] = subjects_dir_backup

    return output


def compute_weighted_mean_surface(in_surfaces):
    """weighted_mean make a weighted average at each node of the surface. The weight are defined by a normal
    distribution (centered on the mid surface)

    Args:
        (list of strings) in_surfaces : List of path to the data projected on the 7 surfaces (35 to 65 % of thickness)
            at each nodes)

    Returns:
        (string) Path to the data averaged
    """
    import os

    import nibabel as nib
    import numpy as np

    # coefficient for normal repartition
    coefficient = [0.1034, 0.1399, 0.1677, 0.1782, 0.1677, 0.1399, 0.1034]

    # sample only to get dimension
    sample = nib.load(in_surfaces[0])
    data_normalized = np.zeros(sample.header.get_data_shape())

    if len(in_surfaces) != 7:
        raise Exception(
            f"There should be 7 surfaces at this point of the pipeline, but found {len(in_surfaces)}, something went wrong..."
        )

    for i in range(len(in_surfaces)):
        current_surf = nib.load(in_surfaces[i])
        data_normalized += current_surf.get_fdata(dtype="float32") * coefficient[i]

    # hemisphere name will always be in our case the first 2 letters of the filename
    hemi = os.path.basename(in_surfaces[0])[0:2]
    # data_normalized = np.atleast_3d(data_normalized)
    hemi_projection = nib.MGHImage(
        data_normalized, affine=sample.affine, header=sample.header
    )
    out_surface = "./" + hemi + ".averaged_projection_on_cortical_surface.mgh"
    out_surface = os.path.abspath(out_surface)
    nib.save(hemi_projection, out_surface)

    return out_surface


def project_onto_fsaverage(
    projection, subject_id, caps_dir, session_id, fwhm, is_longitudinal
):
    """fsaverage_projection projects your data into an averaged subject called fsaverage, available in your $SUBJECTS_DIR
    folder. fsaverage and the subject must be in the subject_dir, so a copy of fsaverage is performed if necessary

    Args:
        (string) projection : Path to the projected data onto native subject surface
        (string) subject_id : The subject id (something like sub-ADNI002S4213)
        (string) session_id : The session id ( something like : ses-M012)
        (string) caps_dir   : Path to the CAPS directory
        (float) fwhm        : FWHM of the Gaussian filter used for smoothing on fsaverage surface (not volume !)
        (bool) is_longitudinal : longitudinal pipeline or not

    Returns:
        (string) Path to the data averaged
    """
    import os
    import shutil

    from nipype.interfaces.freesurfer import MRISPreproc
    from pipelines.pet.surface.utils import get_new_subjects_dir

    subjects_dir_backup = os.path.expandvars("$SUBJECTS_DIR")

    root_env, freesurfer_id = get_new_subjects_dir(
        is_longitudinal, caps_dir, subject_id, session_id
    )

    os.environ["SUBJECTS_DIR"] = root_env

    # copy fsaverage folder next to : subject_id + '_' + session_id
    # for the mris_preproc command to properly find src and target
    fsaverage_has_been_copied = False
    if not os.path.exists(
        os.path.join(os.path.expandvars("$SUBJECTS_DIR"), "fsaverage")
    ):
        shutil.copytree(
            os.path.join(subjects_dir_backup, "fsaverage"),
            os.path.join(os.path.expandvars("$SUBJECTS_DIR"), "fsaverage"),
        )
        fsaverage_has_been_copied = True

    # also copy the mgh file in the surf folder (needed by MRISPreproc
    projection_in_surf_folder = os.path.join(
        os.path.expandvars("$SUBJECTS_DIR"),
        freesurfer_id,
        "surf",
        os.path.basename(projection),
    )

    if not os.path.exists(projection_in_surf_folder):
        shutil.copy(projection, projection_in_surf_folder)

    hemi = os.path.basename(projection)[0:2]
    out_fsaverage = os.path.abspath(
        "./fsaverage_fwhm-" + str(fwhm) + "_" + os.path.basename(projection)
    )

    # Use standalone node
    fsproj = MRISPreproc()
    fsproj.inputs.target = "fsaverage"
    fsproj.inputs.subjects = [freesurfer_id]
    fsproj.inputs.fwhm = fwhm
    fsproj.inputs.hemi = hemi
    fsproj.inputs.surf_measure = os.path.basename(projection)[3:]
    fsproj.inputs.out_file = out_fsaverage
    fsproj.run()

    # remove projection file from surf folder
    os.remove(projection_in_surf_folder)

    # remove fsaverage if it has been copied
    if fsaverage_has_been_copied:
        shutil.rmtree(os.path.join(os.path.expandvars("$SUBJECTS_DIR"), "fsaverage"))

    # put back original subjects_dir env
    os.environ["SUBJECTS_DIR"] = subjects_dir_backup
    return out_fsaverage


def get_mid_surface(in_surfaces):
    """get_mid_surface gives the mid surface when dealing with the 7 different surfaces

    Args:
        (list of strings) in_surfaces : List of path to the 7 different surfaces generated by mris_expand

    Returns:
        (string) Path to the mid surface
    """
    return in_surfaces[3]


def reformat_surfname(hemi, left_surface, right_surface):
    if hemi == "lh":
        return left_surface
    if hemi == "rh":
        return right_surface
    raise ValueError(
        f"First input of this reformat_surfname function must be either lh or rh. Here it is : {hemi}"
    )


def produce_tsv(pet, atlas_files):
    """produce_tsv computes the average of PET signal based on annot files from Freesurfer. Those files describes the
    brain according to known atlases.

        Args:
            (string) pet      : list of path to the PET projection (must be a MGH file) [left_hemisphere, right_hemisphere]
            (string) atlas_files  : Dictionary containing path to lh and rh annotation files for any number of atlases.

        Returns:
            (string) tsv  : path to the tsv containing average PET values
    """
    import os

    import nibabel as nib
    import numpy as np
    import pandas as pds

    # Extract data from projected PET data
    lh_pet_mgh = np.squeeze(nib.load(pet[0]).get_fdata(dtype="float32"))
    rh_pet_mgh = np.squeeze(nib.load(pet[1]).get_fdata(dtype="float32"))

    filename_tsv = []
    for atlas in atlas_files:
        annot_atlas_left = nib.freesurfer.io.read_annot(
            atlas_files[atlas]["lh"], orig_ids=False
        )
        annot_atlas_left[0][annot_atlas_left[0] == -1] = 0
        annot_atlas_right = nib.freesurfer.io.read_annot(
            atlas_files[atlas]["rh"], orig_ids=False
        )
        annot_atlas_right[0][annot_atlas_right[0] == -1] = 0

        average_region = []
        region_names = []
        for r in range(len(annot_atlas_left[2])):
            # cprint(annot_atlas_left[2][r])
            region_names.append(annot_atlas_left[2][r].astype(str) + "_lh")
            region_names.append(annot_atlas_left[2][r].astype(str) + "_rh")

            mask_left = annot_atlas_left[0] == r
            mask_left = np.uint(mask_left)

            masked_data_left = mask_left * lh_pet_mgh
            if np.sum(mask_left) == 0:
                average_region.append(np.nan)
            else:
                average_region.append(np.sum(masked_data_left) / np.sum(mask_left))

            mask_right = annot_atlas_right[0] == r
            mask_right = np.uint(mask_right)
            masked_data_right = mask_right * rh_pet_mgh
            if np.sum(mask_right) == 0:
                average_region.append(np.nan)
            else:
                average_region.append(np.sum(masked_data_right) / np.sum(mask_right))

        final_tsv = pds.DataFrame(
            {
                "index": range(len(region_names)),
                "label_name": region_names,
                "mean_scalar": list(average_region),
            }
        )
        filename_atlas_tsv = "./" + atlas + ".tsv"
        filename_tsv.append(filename_atlas_tsv)
        final_tsv.to_csv(
            filename_atlas_tsv,
            sep="\t",
            index=False,
            columns=["index", "label_name", "mean_scalar"],
        )
    return os.path.abspath(filename_tsv[0]), os.path.abspath(filename_tsv[1])


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
