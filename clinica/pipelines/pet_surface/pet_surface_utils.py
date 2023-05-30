def get_new_subjects_dir(is_longitudinal, caps_dir, subject_id, session_id):
    """Extract SUBJECT_DIR.

    Extract path to FreeSurfer segmentation in CAPS folder and FreeSurfer ID
    (e.g. sub-CLNC01_ses-M000.long.sub-CLNC01_long-M000M018 or sub-CLNC01_ses-M000).
    """
    import os

    from clinica.utils.exceptions import ClinicaCAPSError

    if is_longitudinal:
        root = os.path.join(caps_dir, "subjects", subject_id, session_id, "t1")
        long_folds = [f for f in os.listdir(root) if f.startswith("long-")]
        if len(long_folds) > 1:
            raise ClinicaCAPSError(
                f"[Error] Folder {root} contains {len(long_folds)} folders labeled long-*. Only 1 can exist."
            )
        elif len(long_folds) == 0:
            raise ClinicaCAPSError(
                f"[Error] Folder {root} does not contains a folder labeled long-*. Have you run t1-freesurfer-longitudinal?"
            )
        else:
            root_env = os.path.join(root, long_folds[0], "freesurfer_longitudinal")
            sub_id_cmd = f"{subject_id}_{session_id}.long.{subject_id}_{long_folds[0]}"
    else:
        root_env = os.path.join(
            caps_dir,
            "subjects",
            subject_id,
            session_id,
            "t1",
            "freesurfer_cross_sectional",
        )
        sub_id_cmd = subject_id + "_" + session_id

    return root_env, sub_id_cmd


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

    import clinica.pipelines.pet_surface.pet_surface_utils as utils

    # Old subject_dir is saved for later
    subjects_dir_backup = os.path.expandvars("$SUBJECTS_DIR")

    root_env, freesurfer_id = utils.get_new_subjects_dir(
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

    from clinica.utils.stream import cprint

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

    reg = list(convert_lut.REGION)

    # Check that each label of original volume (old_label) has a matching transformation in the csv file
    for i in range(old_labels.size):
        index = numpy.argwhere(src_val == old_labels[i])
        # Size 0 means no occurrence found
        if index.size == 0:
            raise Exception(
                f"Could not find label {old_labels[i]} on conversion table. Add it manually in CSV file to correct error"
            )
        # else:
        #     cprint(str(old_labels[i]) + " has been found on conversion table")

    # Instantiation of final volume, with same dtype as original volume
    new_volume = numpy.zeros(volume.shape, dtype=volume.dtype)
    # Computing the transformation
    for i in range(len(src)):
        new_volume[volume == src_val[i]] = dst_val[i]
        # cprint("Region " + reg[i] + " transformed")
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
        # cprint("Label " + str(new_labels[i]) + " created in " + current_path)
        allsum = allsum + region_volume

    # The sum of a voxel location across the fourth dimension should be 1
    sum_voxel_mean = float(sum(sum(sum(allsum)))) / allsum.size
    if not isclose(1.0, sum_voxel_mean):
        raise Exception(
            f"Problem during parcellation: the mean sum of a voxel across 4th dimension is {sum_voxel_mean} instead of 1.0"
        )
    # The list of files is returned
    return list_of_regions


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
    # TODO carreful here...
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


def get_wf(
    subject_id,
    session_id,
    pvc_psf_tsv,
    caps_dir,
    pet,
    orig_nu,
    white_surface_left,
    white_surface_right,
    working_directory_subjects,
    acq_label,
    csv_segmentation,
    suvr_reference_region,
    matscript_folder_inverse_deformation,
    destrieux_left,
    destrieux_right,
    desikan_left,
    desikan_right,
    spm_standalone_is_available,
    is_longitudinal,
):
    """get_wf create a full workflow for only one subject, and then executes it

    Args:
        subject_id (string): The subject ID
        session_id (string): The session ID
        pvc_psf_tsv (string): Path the TSV file containing information on the point spread function (PSF)
        caps_dir (string): Path to the CAPS directory
        pet (string): Path to the PET image in the bids directory
        orig_nu (string): Path to the orig_nu file (must be in the CAPS directory, in mri)
        white_surface_left (string): Path to the left white surface in native space of subject
        white_surface_right (string): Path to the right white surface in native space of subject
        working_directory_subjects (string):
        acq_label (string):
        csv_segmentation (string): Path to the CSV for the segmentation (problems encountered while using __file__)
        suvr_reference_region (string): Label of the SUVR reference region
        matscript_folder_inverse_deformation (string): Path to the current folder containing the matlab script used to call the spm function for the inverse deformation
        destrieux_left (string):
        destrieux_right (string):
        desikan_left (string):
        desikan_right (string):
        spm_standalone_is_available (string):
        is_longitudinal (string):

    Returns:
        Void
    """
    import os

    import nipype.interfaces.io as nio
    import nipype.interfaces.utility as niu
    import nipype.pipeline.engine as pe
    from nipype.interfaces.freesurfer import ApplyVolTransform, MRIConvert, Tkregister2
    from nipype.interfaces.fsl import Merge
    from nipype.interfaces.petpvc import PETPVC
    from nipype.interfaces.spm import Coregister, Normalize12

    import clinica.pipelines.pet_surface.pet_surface_utils as utils
    from clinica.utils.filemanip import get_subject_id, load_volume, unzip_nii
    from clinica.utils.pet import get_suvr_mask, read_psf_information
    from clinica.utils.spm import get_tpm
    from clinica.utils.ux import print_begin_image

    image_id = get_subject_id(pet)
    try:
        load_volume(pet)
    except ValueError as e:
        raise ValueError(
            f"Clinica could not load volumes for {image_id.replace('_', ' | ')}. {str(e)}"
        )
    print_begin_image(image_id)

    # Creation of workflow
    # 1 Creation of node
    unzip_pet = pe.Node(
        niu.Function(
            input_names=["in_file"], output_names=["out_file"], function=unzip_nii
        ),
        name="unzip_pet",
    )

    unzip_orig_nu = unzip_pet.clone(name="unzip_orig_nu")

    unzip_mask = unzip_pet.clone(name="unzip_mask")
    unzip_mask.inputs.in_file = get_suvr_mask(suvr_reference_region)

    coreg = pe.Node(Coregister(), name="coreg")

    convert_mgh = pe.Node(MRIConvert(), name="convert_mgh")

    removenan = pe.Node(
        niu.Function(
            input_names=["volname"],
            output_names=["vol_wo_nan"],
            function=utils.remove_nan,
        ),
        name="removenan",
    )

    gtmsegmentation = pe.Node(
        niu.Function(
            input_names=["caps_dir", "subject_id", "session_id", "is_longitudinal"],
            output_names=["gtmseg_file"],
            function=utils.perform_gtmseg,
        ),
        name="gtmseg",
    )
    gtmsegmentation.inputs.caps_dir = caps_dir
    gtmsegmentation.inputs.subject_id = subject_id
    gtmsegmentation.inputs.session_id = session_id
    gtmsegmentation.inputs.is_longitudinal = is_longitudinal

    tkregister = pe.Node(Tkregister2(reg_header=True), name="tkreg")

    convert_gtmseg = convert_mgh.clone(name="convert_gtmseg")

    labelconversion = pe.Node(
        niu.Function(
            input_names=["gtmsegfile", "csv"],
            output_names=["list_of_regions"],
            function=utils.make_label_conversion,
        ),
        name="conversion_of_labels",
    )

    labelconversion.inputs.csv = csv_segmentation

    if not os.path.exists(labelconversion.inputs.csv):
        raise Exception("CSV file : " + labelconversion.inputs.csv + " does not exist.")

    merge_volume = pe.Node(
        Merge(output_type="NIFTI_GZ", dimension="t"), name="merge_volume"
    )

    vol2vol = pe.Node(
        ApplyVolTransform(reg_header=True, interp="trilin"), name="vol2vol"
    )

    vol2vol_mask = pe.Node(
        ApplyVolTransform(reg_header=True, interp="nearest"), name="vol2vol_mask"
    )

    normalize12 = pe.Node(
        Normalize12(
            tpm=get_tpm(),
            affine_regularization_type="mni",
            jobtype="est",
            bias_fwhm=60,
            bias_regularization=0.0001,
            warping_regularization=[0, 0.001, 0.5, 0.05, 0.2],
        ),
        name="normalize_to_MNI",
    )

    if spm_standalone_is_available:
        fun_apply_inverse_deformation = (
            utils.runApplyInverseDeformationField_SPM_standalone
        )
    else:
        fun_apply_inverse_deformation = utils.runApplyInverseDeformationField
    apply_inverse_deformation = pe.Node(
        niu.Function(
            input_names=["target", "deformation_field", "img", "matscript_folder"],
            output_names=["freesurfer_space_eroded_mask"],
            function=fun_apply_inverse_deformation,
        ),
        name="applyInverseDeformation",
    )
    apply_inverse_deformation.inputs.matscript_folder = (
        matscript_folder_inverse_deformation
    )

    pons_normalization = pe.Node(
        niu.Function(
            input_names=["pet_path", "mask"],
            output_names=["suvr"],
            function=utils.suvr_normalization,
        ),
        name="pons_normalization",
    )
    pons_normalization.inputs.pet_tracer = acq_label

    # read_psf_information expects a list of subjects/sessions and returns a list of PSF
    psf_info = read_psf_information(pvc_psf_tsv, [subject_id], [session_id], acq_label)[
        0
    ]

    pvc = pe.Node(PETPVC(pvc="IY"), name="petpvc")
    pvc.inputs.fwhm_x = psf_info[0]
    pvc.inputs.fwhm_y = psf_info[1]
    pvc.inputs.fwhm_z = psf_info[2]

    reformat_surface_name = pe.Node(
        niu.Function(
            input_names=["hemi", "left_surface", "right_surface"],
            output_names=["out"],
            function=utils.reformat_surfname,
        ),
        name="reformat_surface_name",
    )
    reformat_surface_name.inputs.left_surface = white_surface_left
    reformat_surface_name.inputs.right_surface = white_surface_right
    reformat_surface_name.iterables = ("hemi", ["lh", "rh"])

    mris_exp = pe.Node(
        niu.Function(
            input_names=["in_surface"],
            output_names=["out_surface"],
            function=utils.mris_expand,
        ),
        name="mris_expand_white",
    )

    surf_conversion = pe.MapNode(
        niu.Function(
            input_names=[
                "in_surface",
                "reg_file",
                "gtmsegfile",
                "subject_id",
                "session_id",
                "caps_dir",
                "is_longitudinal",
            ],
            output_names=["tval"],
            function=utils.surf2surf,
        ),
        name="surf_conversion",
        iterfield=["in_surface"],
    )
    surf_conversion.inputs.subject_id = subject_id
    surf_conversion.inputs.session_id = session_id
    surf_conversion.inputs.caps_dir = caps_dir
    surf_conversion.inputs.is_longitudinal = is_longitudinal

    vol_on_surf = pe.MapNode(
        niu.Function(
            input_names=[
                "volume",
                "surface",
                "subject_id",
                "session_id",
                "caps_dir",
                "gtmsegfile",
                "is_longitudinal",
            ],
            output_names=["output"],
            function=utils.vol2surf,
        ),
        name="vol_on_surf",
        iterfield=["surface"],
    )
    vol_on_surf.inputs.subject_id = subject_id
    vol_on_surf.inputs.session_id = session_id
    vol_on_surf.inputs.caps_dir = caps_dir
    vol_on_surf.inputs.is_longitudinal = is_longitudinal

    normal_average = pe.Node(
        niu.Function(
            input_names=["in_surfaces"],
            output_names=["out_surface"],
            function=utils.weighted_mean,
        ),
        name="normal_average",
    )

    project_on_fsaverage = pe.Node(
        niu.Function(
            input_names=[
                "projection",
                "subject_id",
                "caps_dir",
                "session_id",
                "fwhm",
                "is_longitudinal",
            ],
            output_names=["out_fsaverage"],
            function=utils.fsaverage_projection,
        ),
        name="project_on_fsaverage",
    )
    project_on_fsaverage.iterables = ("fwhm", [0, 5, 10, 15, 20, 25])
    project_on_fsaverage.inputs.subject_id = subject_id
    project_on_fsaverage.inputs.session_id = session_id
    project_on_fsaverage.inputs.caps_dir = caps_dir
    project_on_fsaverage.inputs.is_longitudinal = is_longitudinal

    extract_mid_surface = pe.Node(
        niu.Function(
            input_names=["in_surfaces"],
            output_names=["mid_surface"],
            function=utils.get_mid_surface,
        ),
        name="extract_mid_surface",
    )

    surface_atlas = {
        "destrieux": {"lh": destrieux_left, "rh": destrieux_right},
        "desikan": {"lh": desikan_left, "rh": desikan_right},
    }

    gather_pet_projection = pe.JoinNode(
        niu.IdentityInterface(fields=["pet_projection_lh_rh"]),
        name="gather_pet_projection_hemisphere",
        joinsource="reformat_surface_name",
        joinfield=["pet_projection_lh_rh"],
    )

    atlas_tsv = pe.Node(
        niu.Function(
            input_names=["pet", "atlas_files"],
            output_names=["destrieux_tsv", "desikan_tsv"],
            function=utils.produce_tsv,
        ),
        name="atlas_tsv",
    )
    atlas_tsv.inputs.atlas_files = surface_atlas

    # 2 creation of workflow : working dir, inputnode, outputnode and datasink
    name_workflow = subject_id.replace("-", "_") + "_" + session_id.replace("-", "_")
    if is_longitudinal:
        name_workflow += "_long"

    wf = pe.Workflow(name=name_workflow)
    wf.base_dir = working_directory_subjects

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "orig_nu",
                "pet",
                "psf",
                "white_surface_left",
                "white_surface_right",
            ]
        ),
        name="inputnode",
        mandatory_inputs=True,
    )

    inputnode.inputs.orig_nu = orig_nu
    inputnode.inputs.pet = pet
    inputnode.inputs.psf = pvc_psf_tsv
    inputnode.inputs.white_surface_right = white_surface_right
    inputnode.inputs.white_surface_left = white_surface_left

    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "mid_surf",
                "projection_native_subject",
                "projection_fsaverage_smoothed",
                "destrieux_tsv",
                "desikan_tsv",
            ]
        ),
        name="outputnode",
        mandatory_inputs=True,
    )

    datasink = pe.Node(nio.DataSink(), name="sinker")

    def get_output_dir(is_longitudinal, caps_dir, subject_id, session_id):
        import os

        from clinica.utils.exceptions import ClinicaCAPSError

        if is_longitudinal:
            root = os.path.join(caps_dir, "subjects", subject_id, session_id, "t1")
            long_folds = [f for f in os.listdir(root) if f.startswith("long-")]
            if len(long_folds) > 1:
                raise ClinicaCAPSError(
                    f"[Error] Folder {root} contains {len(long_folds)} folders labeled long-*. Only 1 can exist"
                )
            elif len(long_folds) == 0:
                raise ClinicaCAPSError(
                    f"[Error] Folder {root} does not contains a folder labeled long-*. Have you run t1-freesurfer-longitudinal?"
                )
            else:
                output_dir = os.path.join(
                    caps_dir,
                    "subjects",
                    subject_id,
                    session_id,
                    "pet",
                    long_folds[0],
                    "surface_longitudinal",
                )
        else:
            output_dir = os.path.join(
                caps_dir, "subjects", subject_id, session_id, "pet", "surface"
            )

        return output_dir

    datasink.inputs.base_directory = get_output_dir(
        is_longitudinal, caps_dir, subject_id, session_id
    )
    datasink.inputs.parameterization = True
    cross_sectional_regexp_substitutions = [
        # Mid surface
        (
            r"(.*(sub-.*)\/(ses-.*)\/pet\/surface)\/midsurface\/.*_hemi_([a-z]+)(.*)$",
            r"\1/\2_\3_hemi-\4_midcorticalsurface",
        ),
        # Projection in native space
        (
            r"(.*(sub-.*)\/(ses-.*)\/pet\/surface)\/projection_native\/.*_hemi_([a-z]+).*",
            r"\1/\2_\3_trc-"
            + acq_label
            + r"_pet_space-native_suvr-"
            + suvr_reference_region
            + r"_pvc-iy_hemi-\4_projection.mgh",
        ),
        # Projection in fsaverage
        (
            r"(.*(sub-.*)\/(ses-.*)\/pet\/surface)\/projection_fsaverage\/.*_hemi_([a-z]+).*_fwhm_([0-9]+).*",
            r"\1/\2_\3_trc-"
            + acq_label
            + r"_pet_space-fsaverage_suvr-"
            + suvr_reference_region
            + r"_pvc-iy_hemi-\4_fwhm-\5_projection.mgh",
        ),
        # TSV file for Destrieux atlas
        (
            r"(.*(sub-.*)\/(ses-.*)\/pet\/surface)\/destrieux_tsv\/destrieux.tsv",
            r"\1/atlas_statistics/\2_\3_trc-"
            + acq_label
            + "_pet_space-destrieux_pvc-iy_suvr-"
            + suvr_reference_region
            + "_statistics.tsv",
        ),
        # TSV file for Desikan atlas
        (
            r"(.*(sub-.*)\/(ses-.*)\/pet\/surface)\/desikan_tsv\/desikan.tsv",
            r"\1/atlas_statistics/\2_\3_trc-"
            + acq_label
            + "_pet_space-desikan_pvc-iy_suvr-"
            + suvr_reference_region
            + "_statistics.tsv",
        ),
    ]
    longitudinal_regexp_substitutions = [
        # Mid surface
        (
            r"(.*(sub-.*)\/(ses-.*)\/pet\/(long-.*)\/surface_longitudinal)\/midsurface\/.*_hemi_([a-z]+)(.*)$",
            r"\1/\2_\3_\4_hemi-\5_midcorticalsurface",
        ),
        # Projection in native space
        (
            r"(.*(sub-.*)\/(ses-.*)\/pet\/(long-.*)\/surface_longitudinal)\/projection_native\/.*_hemi_([a-z]+).*",
            r"\1/\2_\3_\4_trc-"
            + acq_label
            + r"_pet_space-native_suvr-"
            + suvr_reference_region
            + r"_pvc-iy_hemi-\5_projection.mgh",
        ),
        # Projection in fsaverage
        (
            r"(.*(sub-.*)\/(ses-.*)\/pet\/(long-.*)\/surface_longitudinal)\/projection_fsaverage\/.*_hemi_([a-z]+).*_fwhm_([0-9]+).*",
            r"\1/\2_\3_\4_trc-"
            + acq_label
            + r"_pet_space-fsaverage_suvr-"
            + suvr_reference_region
            + r"_pvc-iy_hemi-\5_fwhm-\6_projection.mgh",
        ),
        # TSV file for Destrieux atlas
        (
            r"(.*(sub-.*)\/(ses-.*)\/pet\/(long-.*)\/surface_longitudinal)\/destrieux_tsv\/destrieux.tsv",
            r"\1/atlas_statistics/\2_\3_\4_trc-"
            + acq_label
            + "_pet_space-destrieux_pvc-iy_suvr-"
            + suvr_reference_region
            + "_statistics.tsv",
        ),
        # TSV file for Desikan atlas
        (
            r"(.*(sub-.*)\/(ses-.*)\/pet\/(long-.*)\/surface_longitudinal)\/desikan_tsv\/desikan.tsv",
            r"\1/atlas_statistics/\2_\3_\4_trc-"
            + acq_label
            + "_pet_space-desikan_pvc-iy_suvr-"
            + suvr_reference_region
            + "_statistics.tsv",
        ),
    ]
    if is_longitudinal:
        datasink.inputs.regexp_substitutions = longitudinal_regexp_substitutions
    else:
        datasink.inputs.regexp_substitutions = cross_sectional_regexp_substitutions

    # 3 Connecting the nodes

    # TODO(@arnaud.marcoux): Add titles for sections of connections.
    #   Could be useful to add sections title to group similar connections
    #   together.

    # fmt: off
    wf.connect(
        [
            (inputnode, unzip_pet, [("pet", "in_file")]),
            (unzip_pet, coreg, [("out_file", "source")]),
            (inputnode, convert_mgh, [("orig_nu", "in_file")]),
            (convert_mgh, unzip_orig_nu, [("out_file", "in_file")]),
            (unzip_orig_nu, coreg, [("out_file", "target")]),
            (coreg, removenan, [("coregistered_source", "volname")]),
            (removenan, vol2vol, [("vol_wo_nan", "source_file")]),
            (inputnode, tkregister, [("orig_nu", "target_image")]),
            (unzip_orig_nu, normalize12, [("out_file", "image_to_align")]),
            (unzip_mask, apply_inverse_deformation, [("out_file", "img")]),
            (normalize12, apply_inverse_deformation, [("deformation_field", "deformation_field")]),
            (unzip_orig_nu, apply_inverse_deformation, [("out_file", "target")]),
            (apply_inverse_deformation, vol2vol_mask, [("freesurfer_space_eroded_mask", "source_file")]),
            (gtmsegmentation, vol2vol_mask, [("gtmseg_file", "target_file")]),
            (gtmsegmentation, tkregister, [("gtmseg_file", "moving_image")]),
            (gtmsegmentation, convert_gtmseg, [("gtmseg_file", "in_file")]),
            (gtmsegmentation, vol2vol, [("gtmseg_file", "target_file")]),
            (vol2vol, pons_normalization, [("transformed_file", "pet_path")]),
            (vol2vol_mask, pons_normalization, [("transformed_file", "mask")]),
            (convert_gtmseg, labelconversion, [("out_file", "gtmsegfile")]),
            (labelconversion, merge_volume, [("list_of_regions", "in_files")]),
            (merge_volume, pvc, [("merged_file", "mask_file")]),
            (pons_normalization, pvc, [("suvr", "in_file")]),
            (reformat_surface_name, mris_exp, [("out", "in_surface")]),
            (mris_exp, extract_mid_surface, [("out_surface", "in_surfaces")]),
            (mris_exp, surf_conversion, [("out_surface", "in_surface")]),
            (tkregister, surf_conversion, [("reg_file", "reg_file")]),
            (gtmsegmentation, surf_conversion, [("gtmseg_file", "gtmsegfile")]),
            (pvc, vol_on_surf, [("out_file", "volume")]),
            (surf_conversion, vol_on_surf, [("tval", "surface")]),
            (gtmsegmentation, vol_on_surf, [("gtmseg_file", "gtmsegfile")]),
            (vol_on_surf, normal_average, [("output", "in_surfaces")]),
            (normal_average, project_on_fsaverage, [("out_surface", "projection")]),
            (normal_average, gather_pet_projection, [("out_surface", "pet_projection_lh_rh")]),
            (gather_pet_projection, atlas_tsv, [("pet_projection_lh_rh", "pet")]),
            (atlas_tsv, outputnode, [("destrieux_tsv", "destrieux_tsv")]),
            (atlas_tsv, outputnode, [("desikan_tsv", "desikan_tsv")]),
            (project_on_fsaverage, outputnode, [("out_fsaverage", "projection_fsaverage_smoothed")]),
            (extract_mid_surface, outputnode, [("mid_surface", "mid_surf")]),
            (normal_average, outputnode, [("out_surface", "projection_native_subject")]),
            (outputnode, datasink, [("projection_fsaverage_smoothed", "projection_fsaverage")]),
            (outputnode, datasink, [("mid_surf", "midsurface")]),
            (outputnode, datasink, [("projection_native_subject", "projection_native")]),
            (outputnode, datasink, [("destrieux_tsv", "destrieux_tsv")]),
            (outputnode, datasink, [("desikan_tsv", "desikan_tsv")]),
        ]
    )
    # wf.write_graph(graph2use='flat')
    wf.run()
    # fmt: on
