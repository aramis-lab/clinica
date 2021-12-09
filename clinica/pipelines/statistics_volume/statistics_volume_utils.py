def get_group_1_and_2(tsv, contrast):
    """
        Based on the TSV file given in parameter, compute indexes of each group
    Args:
        tsv: (str) path to the tsv file containing information on subjects/sessions with all covariates
        contrast: (str) name of a column of the tsv

    Returns:
        first_group_idx: (list of int) list of indexes of first group
        second_group_idx: (list of int) list of indexes of second group
        class_names: (list of str of len 2) list of the class names read in the column contrast of the tsv
    """
    import pandas as pds

    from clinica.utils.exceptions import ClinicaException

    # StatisticsVolume pipeline has been instantiated with tsv_file=tsv,
    # so check of existence and integrity have succeeded
    # No further check are done when trying to read it
    tsv = pds.read_csv(tsv, sep="\t")
    columns = list(tsv.columns)

    # An error is raised if the contrast column is not found in the tsv
    if contrast not in columns:
        raise ClinicaException(contrast + " is not present in " + tsv)

    # list(set(my_list)) gives unique values of my_list
    class_names = list(set(tsv[contrast]))

    # This is a 2-sample t-test: we can only allow 2 classes
    if len(class_names) != 2:
        raise ClinicaException(
            "It must exist only 2 classes in the column "
            + contrast
            + " to perform 2-sample t-tests. Here Clinica found: "
            + str(class_names)
        )
    first_group_idx = [
        i for i, label in enumerate(list(tsv[contrast])) if label == class_names[0]
    ]
    second_group_idx = [
        i for i, label in enumerate(list(tsv[contrast])) if label == class_names[1]
    ]

    return first_group_idx, second_group_idx, class_names


def model_creation(tsv, contrast, idx_group1, idx_group2, file_list, template_file):
    """
        Create the matlab .m file for the instantiation of the 2-sample t-test model in SPM

    Args:
        tsv: (str) path to the tsv file containing information on subjects/sessions with all covariates
        contrast: (str) name of a column of the tsv
        idx_group1: (list of int) list of indexes of first group
        idx_group2: (list of int) list of indexes of second group
        file_list: List of files used in the statistical test. Their order is the same as it appears on the tsv file
        template_file: (str) path to the template file used to generate the .m file

    Returns:
        current_model: (str) path to the matlab files with all the @TEXT replaced with the correct names
        covariates: list of str with the names of covariates
    """
    from numbers import Number
    from os import mkdir, remove
    from os.path import abspath, dirname, isdir, isfile, join
    from shutil import rmtree

    import numpy as np
    import pandas as pds

    import clinica.pipelines.statistics_volume.statistics_volume_utils as utls
    from clinica.utils.exceptions import ClinicaException

    # Get template for model creation
    if not isfile(template_file):
        raise RuntimeError("[Error] " + template_file + " could not be found !")

    # Create current model filename
    current_model = abspath("./current_model_creation.m")

    if isfile(current_model):
        # This should never be reached
        remove(current_model)

    # Read template
    with open(template_file, "r") as file:
        filedata = file.read()

    # Create output folder for SPM (and remove it if it already exists)
    output_folder = abspath(join(dirname(current_model), "..", "2_sample_t_test"))
    if isdir(output_folder):
        rmtree(output_folder)
    mkdir(output_folder)

    # Replace string in matlab file to set the output directory, and all the scans used for group 1 and 2
    filedata = filedata.replace("@OUTPUTDIR", "'" + output_folder + "'")
    filedata = filedata.replace(
        "@SCANS1",
        utls.unravel_list_for_matlab(
            [f for i, f in enumerate(file_list) if i in idx_group1]
        ),
    )
    filedata = filedata.replace(
        "@SCANS2",
        utls.unravel_list_for_matlab(
            [f for i, f in enumerate(file_list) if i in idx_group2]
        ),
    )

    # Write filedata in the output file
    with open(current_model, "w+") as file:
        file.write(filedata)

    # Add our covariates
    tsv = pds.read_csv(tsv, sep="\t")
    columns_stripped = [elem.strip(" ") for elem in list(tsv.columns)]
    if columns_stripped != list(tsv.columns):
        raise ClinicaException(
            "[Error] Check the column of your tsv file "
            + tsv
            + "Whitespace in the column names can cause errors"
        )
    covariates = [
        elem
        for elem in columns_stripped
        if elem not in ["participant_id", "session_id", contrast]
    ]
    for covar_number, covar in enumerate(covariates, start=1):
        current_covar_data = list(tsv[covar])
        if isinstance(current_covar_data[0], str):
            # Transform data
            temp_data = [elem.replace(",", ".") for elem in current_covar_data]
            if all(utls.is_number(elem) for elem in temp_data):
                current_covar_data = [float(elem) for elem in temp_data]
            else:
                # categorical variables (like Male; Female; M, F etc...)
                unique_values = list(np.unique(np.array(current_covar_data)))
                current_covar_data = [
                    unique_values.index(elem) for elem in current_covar_data
                ]

        elif isinstance(current_covar_data[0], Number):
            # Do nothing
            pass

        current_covar_data_group1 = [
            elem for i, elem in enumerate(current_covar_data) if i in idx_group1
        ]
        current_covar_data_group2 = [
            elem for i, elem in enumerate(current_covar_data) if i in idx_group2
        ]
        covar_data_concatenated = current_covar_data_group1 + current_covar_data_group2
        utls.write_covariate_lines(
            current_model, covar_number, covar, covar_data_concatenated
        )

    # Tell matlab to run the script at the end
    with open(current_model, "a") as file:
        file.write("spm_jobman('run', matlabbatch)")
    return current_model, covariates


def is_number(s: str):
    """

    Args:
        s: (str) string to test if it can be converted into float

    Returns:
        True or False
    """
    try:
        # Try the conversion, if it is not possible, error will be raised
        float(s)
        return True
    except ValueError:
        return False


def unravel_list_for_matlab(my_list):
    """
        Unpack a list into a Matlab compliant format to insert in .m files.
    Args:
        my_list: (list) of str

    Returns:
        (str) that join the different str of the list, with '
    """
    result = "','".join(my_list)
    result = "'" + result + "'"
    return result


def write_covariate_lines(m_file_to_write_in, covar_number, covar_name, covar_values):
    """
        Use this function to add covariate lines in the Matlab file m_file_to_write_in for one covariate
    Args:
        m_file_to_write_in: (str) path to the m-file
        covar_number: (int) this is the number of the covariate (must start at 1)
        covar_name: (str) name of the covariate (ex: 'age', 'sex')
        covar_values: (list) of float with the values of the covariates

    Returns:
        nothing
    """
    from os.path import isfile

    assert isfile(m_file_to_write_in), "Could not find file " + m_file_to_write_in
    covar_values_string = [str(elem) for elem in covar_values]
    covar_values_for_script = "[" + " ".join(covar_values_string) + "]"

    m_file = open(m_file_to_write_in, "a")
    m_file.write(
        "matlabbatch{1}.spm.stats.factorial_design.cov("
        + str(covar_number)
        + ").c = "
        + covar_values_for_script
        + ";\n"
    )
    m_file.write(
        "matlabbatch{1}.spm.stats.factorial_design.cov("
        + str(covar_number)
        + ").cname = '"
        + covar_name
        + "';\n"
    )
    m_file.write(
        "matlabbatch{1}.spm.stats.factorial_design.cov("
        + str(covar_number)
        + ").iCFI = 1;\n"
    )
    m_file.write(
        "matlabbatch{1}.spm.stats.factorial_design.cov("
        + str(covar_number)
        + ").iCC = 1;\n"
    )
    m_file.close()


def run_m_script(m_file):
    """
        Runs a matlab m file for SPM, determining automatically if it must be launched with SPM or SPM Standalone
        If launch with spm standalone, the line 'spm_jobman('run', matlabbatch)' must be removed because unnecessary

    Args:
        m_file: (str) path to Matlab m file

    Returns:
        output_mat_file: (str) path to the SPM.mat file needed in SPM analysis
    """
    import platform
    from os import system
    from os.path import abspath, basename, dirname, isfile, join

    from nipype.interfaces.matlab import MatlabCommand, get_matlab_command

    import clinica.pipelines.statistics_volume.statistics_volume_utils as utls
    from clinica.utils.spm import spm_standalone_is_available

    assert isinstance(m_file, str), "[Error] Argument must be a string"
    if not isfile(m_file):
        raise FileNotFoundError("[Error] File " + m_file + "does not exist")
    assert m_file[-2:] == ".m", (
        "[Error] " + m_file + " is not a Matlab file (extension must be .m)"
    )

    # Generate command line to run
    if spm_standalone_is_available():
        utls.delete_last_line(m_file)
        # SPM standalone must be run directly from its root folder
        if platform.system().lower().startswith("darwin"):
            # Mac OS
            cmdline = (
                "cd $SPMSTANDALONE_HOME && ./run_spm12.sh $MCR_HOME batch " + m_file
            )
        elif platform.system().lower().startswith("linux"):
            # Linux OS
            cmdline = "$SPMSTANDALONE_HOME/run_spm12.sh $MCR_HOME batch " + m_file
        else:
            raise SystemError("Clinica only support Mac OS and Linux")
        system(cmdline)
    else:
        MatlabCommand.set_default_matlab_cmd(get_matlab_command())
        matlab = MatlabCommand()
        if platform.system().lower().startswith("linux"):
            matlab.inputs.args = "-nosoftwareopengl"
        matlab.inputs.paths = dirname(m_file)
        matlab.inputs.script = basename(m_file)[:-2]
        matlab.inputs.single_comp_thread = False
        matlab.inputs.logfile = abspath("./matlab_output.log")
        matlab.run()
    output_mat_file = abspath(join(dirname(m_file), "..", "2_sample_t_test", "SPM.mat"))
    if not isfile(output_mat_file):
        raise RuntimeError("Output matrix " + output_mat_file + " was not produced")
    return output_mat_file


def delete_last_line(filename):
    """
    Use this function to remove the call to spm jobman if m file is used with SPM standalone
    Args:
        filename: path to filename

    Returns:
        Nothing
    """
    import os

    with open(filename, "r+", encoding="utf-8") as file:

        # Move the pointer (similar to a cursor in a text editor) to the end of the file
        file.seek(0, os.SEEK_END)

        # This code means the following code skips the very last character in the file -
        # i.e. in the case the last line is null we delete the last line
        # and the penultimate one
        pos = file.tell() - 1

        # Read each character in the file one at a time from the penultimate
        # character going backwards, searching for a newline character
        # If we find a new line, exit the search
        while pos > 0 and file.read(1) != "\n":
            pos -= 1
            file.seek(pos, os.SEEK_SET)

        # So long as we're not at the start of the file, delete all the characters ahead
        # of this position
        if pos > 0:
            file.seek(pos, os.SEEK_SET)
            file.truncate()


def estimate(mat_file, template_file):
    """
        Make a copy of the template file (for estimation) and replace @SPMMAT by the real path to SPM.mat (mat_file)
    Args:
        mat_file: (str) path to the SPM.mat file of the SPM analysis
        template_file: (str) path to the template file for the estimation of the model

    Returns:
        current_model_estimation: (str) path to the template file filled with SPM.mat information, ready to be launched
    """
    from os.path import abspath

    # Read template
    with open(template_file, "r") as file:
        filedata = file.read()
    # Replace by the real path to spm.mat
    filedata = filedata.replace("@SPMMAT", "'" + mat_file + "'")
    current_model_estimation = abspath("./current_model_creation.m")
    with open(current_model_estimation, "w+") as file:
        file.write(filedata)

    return current_model_estimation


def results(mat_file, template_file, method, threshold):
    """
        Make a copy of the template file (for results) and replace @SPMMAT by the real path to SPM.mat (mat_file)
    Args:
        mat_file: (str) path to the SPM.mat file of the SPM analysis
        template_file: (str) path to the template file for getting the results of the model
        method: (str)
        threshold

    Returns:
        current_model_estimation: (str) path to the template file filled with SPM.mat information, ready to be launched
    """
    import time
    from os.path import abspath

    # Read template
    with open(template_file, "r") as file:
        filedata = file.read()
    # Replace by the real path to spm.mat
    filedata = filedata.replace("@SPMMAT", "'" + mat_file + "'")
    filedata = filedata.replace("@CORRECTIONMETHOD", "'" + method + "'")
    filedata = filedata.replace("@THRESH", str(threshold))
    if method != "none":
        # Necessary so that the 2 computations does not overlap
        time.sleep(15)
    current_model_result = abspath("./current_model_result.m")
    with open(current_model_result, "w+") as file:
        file.write(filedata)

    return current_model_result


def contrast(mat_file, template_file, covariates, class_names):
    """
        Make a copy of the template file (for results) and replace @SPMMAT, @COVARNUMBER, @GROUP1, @GROUP2
        by the corresponding variables

    Args:
        mat_file: (str) path to the SPM.mat file of the SPM analysis
        template_file: (str) path to the template file for getting the results of the model
        covariates: (list) of str: list of covariates
        class_names: (list) of str of length 2 that correspond to the 2 classes for the group comparison

    Returns:
        current_model_estimation: (str) path to the template file filled with variables, ready to be launched
    """
    from os.path import abspath

    number_of_covariates = len(covariates)

    # Read template
    with open(template_file, "r") as file:
        filedata = file.read()
    # Replace by the real path to spm.mat
    filedata = filedata.replace("@SPMMAT", "'" + mat_file + "'")
    filedata = filedata.replace("@COVARNUMBER", "0 " * number_of_covariates)
    filedata = filedata.replace("@GROUP1", "'" + class_names[0] + "'")
    filedata = filedata.replace("@GROUP2", "'" + class_names[1] + "'")
    current_model_estimation = abspath("./current_model_contrast.m")
    with open(current_model_estimation, "w+") as file:
        file.write(filedata)

    return current_model_estimation


def read_output(spm_mat, class_names, covariates, group_label, fwhm, measure):
    """
        Once analysis is done, grab all the different filenames and rename them in current directory according to class
        names
    Args:
        spm_mat: (str) path to the SPM.mat file of the SPM analysis
        class_names: (list) of str of length 2 that correspond to the 2 classes for the group comparison
        covariates: (list) of str: list of covariates
        group_label: name of the group label
        fwhm: fwhm in mm used
        measure: measure used

    Returns:
        spmT_0001: (str) path to t maps for the first group comparison
        spmT_0002: (str) path to t maps for the second group comparison
        spm_figures: (list) path to figure files
        variance_of_error: (str) path to variance of error
        resels_per_voxels: (str) path to resels per voxel
        mask: (str) path to mask of included voxels
        regression_coeff: (str list) path to regression coefficients
        contrasts: (str list) path to weighted parameter estimation for the 2 contrasts

    """
    from os import listdir
    from os.path import abspath, dirname, isdir, isfile, join
    from shutil import copyfile

    if not isfile(spm_mat):
        if not isdir(dirname(spm_mat)):
            raise RuntimeError(
                "[Error] output folder " + dirname(spm_mat) + " does not exist."
            )
        else:
            raise RuntimeError("[Error] SPM matrix " + spm_mat + " does not exist.")

    list_files = [f for f in listdir(dirname(spm_mat)) if not f.startswith(".")]

    # Handle figure files
    figures = [
        abspath(join(dirname(spm_mat), f)) for f in list_files if f.endswith("png")
    ]
    if len(figures) < 2:
        raise RuntimeError("[Error] Figures were not generated")
    fig_number = [int(f[-7:-4]) for f in figures]
    spm_figures = [
        abspath("./" + "group-" + group_label + "_report-" + str(i) + ".png")
        for i in fig_number
    ]
    for old_name, new_name in zip(figures, spm_figures):
        copyfile(old_name, new_name)

    # Handle spm t maps
    spm_T = [
        abspath(join(dirname(spm_mat), f)) for f in list_files if f.startswith("spmT")
    ]
    spm_T = sorted(spm_T)
    if len(spm_T) != 2:
        raise RuntimeError("[Error] " + str(len(spm_T)) + " SPM t-map(s) were found")
    if fwhm:
        spmT_0001 = abspath(
            f"group-{group_label}_{class_names[0]}-lt-{class_names[1]}_measure-{measure}_fwhm-{int(fwhm)}_TStatistics.nii"
        )
        spmT_0002 = abspath(
            f"group-{group_label}_{class_names[1]}-lt-{class_names[0]}_measure-{measure}_fwhm-{int(fwhm)}_TStatistics.nii"
        )
    else:
        spmT_0001 = abspath(
            f"group-{group_label}_{class_names[0]}-lt-{class_names[1]}_measure-{measure}_TStatistics.nii"
        )
        spmT_0002 = abspath(
            f"group-{group_label}_{class_names[1]}-lt-{class_names[0]}_measure-{measure}_TStatistics.nii"
        )
    copyfile(join(dirname(spm_mat), "spmT_0001.nii"), spmT_0001)
    copyfile(join(dirname(spm_mat), "spmT_0002.nii"), spmT_0002)

    variance_of_error = abspath("./group-" + group_label + "_VarianceError.nii")
    copyfile(abspath(join(dirname(spm_mat), "ResMS.nii")), variance_of_error)

    resels_per_voxels = abspath("./resels_per_voxel.nii")
    copyfile(abspath(join(dirname(spm_mat), "RPV.nii")), resels_per_voxels)

    mask = abspath("./included_voxel_mask.nii")
    copyfile(abspath(join(dirname(spm_mat), "mask.nii")), mask)

    # Handle beta files
    betas = [
        abspath(join(dirname(spm_mat), f)) for f in list_files if f.startswith("beta_")
    ]
    betas = sorted(betas)
    if len(betas) != 2 + len(covariates):
        raise RuntimeError("[Error] Not enough betas files found in output directory")
    regression_coeff_covar = [abspath("./" + covar + ".nii") for covar in covariates]
    regression_coeff = [
        abspath("./" + class_names[0] + ".nii"),
        abspath("./" + class_names[1] + ".nii"),
    ]
    regression_coeff.extend(regression_coeff_covar)
    # Order is respected:
    for beta, reg_coeff in zip(betas, regression_coeff):
        copyfile(beta, reg_coeff)

    # Handle contrast files
    con_files = [
        abspath(join(dirname(spm_mat), f)) for f in list_files if f.startswith("con_")
    ]
    if len(con_files) != 2:
        raise RuntimeError("There must exists only 2 contrast files !")
    contrasts = [
        abspath(
            f"group-{group_label}_{class_names[0]}-lt-{class_names[1]}_measure-{measure}_contrast.nii"
        ),
        abspath(
            f"group-{group_label}_{class_names[1]}-lt-{class_names[0]}_measure-{measure}_contrast.nii"
        ),
    ]
    for con, contrast in zip(con_files, contrasts):
        copyfile(con, contrast)

    return (
        spmT_0001,
        spmT_0002,
        spm_figures,
        variance_of_error,
        resels_per_voxels,
        mask,
        regression_coeff,
        contrasts,
    )
