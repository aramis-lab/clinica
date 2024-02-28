def runApplyInverseDeformationField_SPM_standalone(
    target, deformation_field, img, matscript_folder
):
    """
    This function does the exact same job as runApplyInverseDeformationField but with SPM standalone. We directly create
    a batch file that SPM standalone can run. This function does not check whether SPM standalone must be used. Previous
    check when building the pipeline ensures that all the env vars exists ($SPMSTANDALONE_HOME and $MCR_HOME)
    """
    import os
    import platform
    import subprocess
    from os.path import abspath, basename, exists, join

    prefix = "subject_space_"

    # Write SPM batch command directly in a script that is readable by SPM standalone
    script_location = abspath("./m_script.m")
    script_file = open(script_location, "w+")
    script_file.write(
        "jobs{1}.spm.util.defs.comp{1}.inv.comp{1}.def = {'"
        + deformation_field
        + "'};\n"
    )
    script_file.write("jobs{1}.spm.util.defs.comp{1}.inv.space = {'" + target + "'};\n")
    script_file.write("jobs{1}.spm.util.defs.out{1}.pull.fnames = {'" + img + "'};\n")
    script_file.write(
        "jobs{1}.spm.util.defs.out{1}.pull.savedir.saveusr = {'"
        + abspath(os.getcwd())
        + "'};\n"
    )
    script_file.write("jobs{1}.spm.util.defs.out{1}.pull.interp = 4;\n")
    script_file.write("jobs{1}.spm.util.defs.out{1}.pull.mask = 1;\n")
    script_file.write("jobs{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];\n")
    script_file.write("jobs{1}.spm.util.defs.out{1}.pull.prefix = '" + prefix + "';\n")
    script_file.close()

    # Generate command line to run
    # SPM standalone must be run directly from its root folder
    if platform.system() == "Darwin":
        # Mac OS
        cmdline = f"cd $SPMSTANDALONE_HOME && ./run_spm12.sh $MCR_HOME batch {script_location}"
    elif platform.system() == "Linux":
        # Linux OS
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


def runApplyInverseDeformationField(target, deformation_field, img, matscript_folder):
    import os
    import sys

    from nipype.interfaces.matlab import MatlabCommand, get_matlab_command

    prefix = "subject_space_"

    MatlabCommand.set_default_matlab_cmd(get_matlab_command())
    matlab = MatlabCommand()
    if sys.platform.startswith("linux"):
        matlab.inputs.args = "-nosoftwareopengl"

    matlab.inputs.paths = matscript_folder
    matlab.inputs.script = """
    applyInverseDeformationField('%s', '%s', '%s', './', '%s')
    """ % (
        target,
        deformation_field,
        img,
        prefix,
    )
    matlab.inputs.single_comp_thread = False
    matlab.inputs.logfile = os.path.join("./", "matlab_output.log")
    matlab.run()

    output_file = os.path.join(os.path.abspath("./"), prefix + os.path.basename(img))
    if not os.path.exists(output_file):
        raise IOError(
            f"Something went wrong, please check {os.path.abspath(matlab.inputs.logfile)} for more information"
        )
    return output_file


def suvr_normalization(pet_path, mask):
    """suvr_normalization is a way of getting suvr from your pet image, based on the segmentation performed by
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
    import shutil
    import subprocess
    import sys

    # bug in mris_expand : you are not allowed to write the surfaces elsewhere than in the surf folder
    # -N is a hidden parameter (not documented) that allows the user to specify the number of surface generated between
    # source and final target surface. Here target is 65% of thickness, with 13 surfaces. Then we only keep the surfaces
    # we are interested in.

    out_file = in_surface + "_exp-"
    cmd = "mris_expand -thickness -N 13 " + in_surface + " 0.65 " + out_file
    # If system is MacOS, this export command must be run just before the mri_vol2surf command to bypass MacOs security
    if sys.platform == "darwin":
        cmd = "export DYLD_LIBRARY_PATH=$FREESURFER_HOME/lib/gcc/lib && " + cmd

    subprocess_mris_expand = subprocess.run(
        cmd,
        shell=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    if subprocess_mris_expand.returncode != 0:
        raise ValueError("mris_expand failed, returned non-zero code")

    # Move surface file to the current directory
    out_filelist = [
        out_file + "007",
        out_file + "008",
        out_file + "009",
        out_file + "010",
        out_file + "011",
        out_file + "012",
        out_file + "013",
    ]
    out_in_node = [os.path.abspath("./" + os.path.basename(x)) for x in out_filelist]
    for i in range(len(out_filelist)):
        shutil.move(out_filelist[i], out_in_node[i])

    # Remove useless surfaces (0%, 5%, 10%, 15%, 20%, 25% and 30% of thickness)
    to_remove = [
        out_file + "000",
        out_file + "001",
        out_file + "002",
        out_file + "003",
        out_file + "004",
        out_file + "005",
        out_file + "006",
    ]
    for i in range(len(to_remove)):
        os.remove(to_remove[i])

    return out_in_node


def surf2surf(
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

    import clinica.pipelines.pet_surface.pet_surface_utils as utils

    # set subjects_dir env. variable for mri_surf2surf to work properly
    subjects_dir_backup = os.path.expandvars("$SUBJECTS_DIR")

    root_env, freesurfer_id = utils.get_new_subjects_dir(
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


def vol2surf(
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

    import clinica.pipelines.pet_surface.pet_surface_utils as utils

    # set subjects_dir env. variable for mri_vol2surf to work properly
    subjects_dir_backup = os.path.expandvars("$SUBJECTS_DIR")

    root_env, freesurfer_id = utils.get_new_subjects_dir(
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


def weighted_mean(in_surfaces):
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


def fsaverage_projection(
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

    import clinica.pipelines.pet_surface.pet_surface_utils as utils

    subjects_dir_backup = os.path.expandvars("$SUBJECTS_DIR")

    root_env, freesurfer_id = utils.get_new_subjects_dir(
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
    res = None
    if hemi == "lh":
        res = left_surface
    elif hemi == "rh":
        res = right_surface
    else:
        raise Exception(
            "First input of this reformat_surfname function must be either lh or rh"
        )
    return res


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

    from clinica.utils.stream import cprint

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
