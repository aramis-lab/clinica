import typing as ty
from pathlib import Path


def get_group_1_and_2(tsv: str, contrast: str) -> ty.Tuple[list]:
    """Based on the TSV file given in parameter, compute indexes of each group

    Parameters
    ----------
    tsv: str
        Path to the tsv file containing information on subjects/sessions with all covariates
    contrast: str
        Name of a column of the tsv

    Returns
    -------
    first_group_idx: list of int
        list of indexes of first group
    second_group_idx: list of int
        list of indexes of second group
    class_names: list of str of len 2
        list of the class names read in the column contrast of the tsv
    """
    import pandas as pds

    from clinica.utils.exceptions import ClinicaException

    # StatisticsVolume pipeline has been instantiated with tsv_file=tsv,
    # so check of existence and integrity have succeeded
    # No further check are done when trying to read it
    tsv = pds.read_csv(tsv, sep="\t")
    columns = list(tsv.columns)

    if contrast not in columns:
        raise ClinicaException(contrast + " is not present in " + tsv)

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


def create_spm_output_folder(current_model: str) -> str:
    """Creates the spm output folder.

    If it already exists, it is deleted and created again.

    Parameters
    ----------
    current_model: str
        Path to the matlab files with all the @TEXT replaced with the correct names

    Returns
    -------
    output_folder: str
        Path to the spm output_folder
    """
    from os import mkdir
    from os.path import abspath, dirname, isdir, join
    from shutil import rmtree

    output_folder = abspath(join(dirname(current_model), "..", "2_sample_t_test"))
    if isdir(output_folder):
        rmtree(output_folder)
    mkdir(output_folder)
    return output_folder


def set_output_and_groups(
    output_folder: str,
    current_model: str,
    file_list: list,
    idx_group1: list,
    idx_group2: list,
    filedata: str,
) -> None:
    """Sets path to output and groups in the .m script

    Parameters
    ----------
    output_folder: str
        Path to the folder where the pipelines output will be generated
    current_model: str
        Path to the matlab files with all the @TEXT replaced with the correct names
    file_list: List
        List of files used in the statistical test. Their order is the same as it appears on the tsv file
    idx_group1: list of int
        List of indexes of first group
    idx_group2: list of int
        List of indexes of second group
    filedata: str
        .m script template which is to be modified and then run
    """
    from clinica.pipelines.statistics_volume.statistics_volume_utils import (
        unravel_list_for_matlab,
    )

    filedata = filedata.replace("@OUTPUTDIR", "'" + output_folder + "'")
    filedata = filedata.replace(
        "@SCANS1",
        unravel_list_for_matlab(
            [f for i, f in enumerate(file_list) if i in idx_group1]
        ),
    )
    filedata = filedata.replace(
        "@SCANS2",
        unravel_list_for_matlab(
            [f for i, f in enumerate(file_list) if i in idx_group2]
        ),
    )

    with open(current_model, "w+") as file:
        file.write(filedata)


def convert_to_numeric(current_covar_data: list) -> list:
    """Convert the string values from the input list to numeric values.

    If the values can be casted to floats, then outputs a list of floats.
    If the values are categorical, then performs encoding of the categories.

    Parameters
    ----------
    current_covar_data: list
    Returns
    -------
    current_covar_data: list
    """
    import numpy as np

    from clinica.pipelines.statistics_volume.statistics_volume_utils import is_number

    temp_data = [elem.replace(",", ".") for elem in current_covar_data]
    if all(is_number(elem) for elem in temp_data):
        return [float(elem) for elem in temp_data]
    # categorical variables (like Male; Female; M, F etc...)
    unique_values = list(np.unique(np.array(current_covar_data)))
    return [unique_values.index(elem) for elem in current_covar_data]


def write_covariates(
    current_covar_data: list,
    idx_group1: list,
    idx_group2: list,
    model: str,
    covar: int,
    covar_number: int,
) -> None:
    """Concatenate covariates for group1 and group2 and write them to disk.

    Parameters
    ----------
    current_covar_data: list
    idx_group1: list of int
        List of indexes of first group
    idx_group2: list of int
        List of indexes of second group
    model: str
        Path to the matlab files with all the @TEXT replaced with the correct names
    covar: int
        Covariance
    covar_number: int
    """
    from clinica.pipelines.statistics_volume.statistics_volume_utils import (
        write_covariate_lines,
    )

    current_covar_data_group1 = [
        elem for i, elem in enumerate(current_covar_data) if i in idx_group1
    ]
    current_covar_data_group2 = [
        elem for i, elem in enumerate(current_covar_data) if i in idx_group2
    ]
    covar_data_concatenated = current_covar_data_group1 + current_covar_data_group2
    write_covariate_lines(model, covar_number, covar, covar_data_concatenated)


def write_matlab_model(
    tsv: str,
    contrast: str,
    idx_group1: list,
    idx_group2: list,
    file_list: list,
    template_file: str,
):
    """Create the matlab .m file for the instantiation of the 2-sample t-test model in SPM

    Parameters
    ----------
    tsv: str
        Path to the tsv file containing information on subjects/sessions with all covariates
    contrast: str
        Name of a column of the tsv
    idx_group1: list of int
        List of indexes of first group
    idx_group2: list of int
        List of indexes of second group
    file_list: List
        List of files used in the statistical test. Their order is the same as it appears on the tsv file
    template_file: str
        Path to the template file used to generate the .m file

    Returns
    -------
    current_model: str
        Path to the matlab files with all the @TEXT replaced with the correct names
    covariates: list of str
        Names of covariates
    """
    from numbers import Number
    from os import remove
    from os.path import abspath, isfile

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
        remove(current_model)

    with open(template_file, "r") as file:
        filedata = file.read()

    output_folder = create_spm_output_folder(current_model)

    # Replace string in matlab file to set the output directory, and all the scans used for group 1 and 2
    set_output_and_groups(
        output_folder, current_model, file_list, idx_group1, idx_group2, filedata
    )

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
            current_covar_data = convert_to_numeric(current_covar_data)
        elif isinstance(current_covar_data[0], Number):
            # Do nothing
            pass
        write_covariates(
            current_covar_data,
            idx_group1,
            idx_group2,
            current_model,
            covar,
            covar_number,
        )

    # Tell matlab to run the script at the end
    with open(current_model, "a") as file:
        file.write("spm_jobman('run', matlabbatch)")
    return current_model, covariates


def is_number(s: str) -> bool:
    """Returns True is the input can be converted to float, False otherwise.

    Parameters
    ----------
    s: str
        String to test if it can be converted into float

    Returns
    -------
    bool
        True or False
    """
    try:
        # Try the conversion, if it is not possible, error will be raised
        float(s)
        return True
    except ValueError:
        return False


def unravel_list_for_matlab(my_list: list) -> str:
    """Unpack a list into a Matlab compliant format to insert in .m files.

    Parameters
    ----------
    my_list: list of str

    Returns
    -------
    str
        Contains the different element of the list joined in one string.
    """
    result = "','".join(my_list)
    result = "'" + result + "'"
    return result


def write_covariate_lines(
    m_file_to_write_in: str, covar_number: int, covar_name: str, covar_values: list
) -> None:
    """Add one line in the Matlab file `m_file_to_write_in` for each covariate.

    Parameters
    ----------
    m_file_to_write_in: str
        Path to the m-file
    covar_number: int
        This is the number of the covariate (must start at 1)
    covar_name: str
        Name of the covariate (ex: 'age', 'sex')
    covar_values: list of float
        Values of the covariates
    """
    from os.path import isfile

    if not isfile(m_file_to_write_in):
        raise FileNotFoundError(f"Could not find file {m_file_to_write_in}")
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


def run_m_script(m_file: str) -> str:
    """Runs a Matlab file for SPM.

    Determines automatically if the script should be launched with SPM or SPM standalone.
    If launch with spm standalone, the line 'spm_jobman('run', matlabbatch)' must be removed because unnecessary

    Parameters
    ----------
    m_file: str
        Path to Matlab m file

    Returns
    -------
    output_mat_file: str
        Path to the SPM.mat file needed in SPM analysis
    """

    from os.path import abspath, dirname, isfile, join
    from pathlib import Path

    from clinica.utils.spm import spm_standalone_is_available

    m_file = Path(m_file)
    if not m_file.exists():
        raise FileNotFoundError("[Error] File " + m_file + "does not exist")
    if m_file.suffix != ".m":
        raise ValueError(
            f"[Error] {m_file} is not a Matlab file (extension must be .m)"
        )

    # Generate command line to run
    if spm_standalone_is_available():
        run_matlab_script_with_spm_standalone(m_file)
    else:
        run_matlab_script_with_matlab(m_file)
    output_mat_file = abspath(join(dirname(m_file), "..", "2_sample_t_test", "SPM.mat"))
    if not isfile(output_mat_file):
        raise RuntimeError("Output matrix " + output_mat_file + " was not produced")
    return output_mat_file


def run_matlab_script_with_matlab(m_file: Path) -> None:
    """Runs a matlab script using matlab

    Parameters
    ----------
    m_file: str
        Path to the script
    """
    import platform
    from os.path import abspath

    from nipype.interfaces.matlab import MatlabCommand, get_matlab_command

    MatlabCommand.set_default_matlab_cmd(get_matlab_command())
    matlab = MatlabCommand()
    if platform.system().lower().startswith("linux"):
        matlab.inputs.args = "-nosoftwareopengl"
    matlab.inputs.paths = m_file.parent
    matlab.inputs.script = m_file.stem
    matlab.inputs.single_comp_thread = False
    matlab.inputs.logfile = abspath("./matlab_output.log")
    matlab.run()


def run_matlab_script_with_spm_standalone(m_file: str) -> None:
    """Runs a matlab script spm_standalone

    Parameters
    ----------
    m_file: str
        Path to the script
    """
    import platform
    from os import system

    import clinica.pipelines.statistics_volume.statistics_volume_utils as utls

    utls.delete_last_line(m_file)
    # SPM standalone must be run directly from its root folder
    if platform.system().lower().startswith("darwin"):
        # Mac OS
        cmdline = "cd $SPMSTANDALONE_HOME && ./run_spm12.sh $MCR_HOME batch " + m_file
    elif platform.system().lower().startswith("linux"):
        # Linux OS
        cmdline = "$SPMSTANDALONE_HOME/run_spm12.sh $MCR_HOME batch " + m_file
    else:
        raise SystemError("Clinica only support Mac OS and Linux")
    system(cmdline)


def delete_last_line(filename: str) -> None:
    """Use this function to remove the call to spm jobman if m file is used with SPM standalone

    Parameters
    ----------
    filename: str
        Path to filename
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


def clean_template_file(mat_file: str, template_file: str) -> str:
    """Make a copy of the template file (for estimation) and replace @SPMMAT by the real path to SPM.mat (mat_file)

    Parameters
    ----------
    mat_file: str
        Path to the SPM.mat file of the SPM analysis
    template_file: str
        Path to the template file for the estimation of the model

    Returns
    -------
    current_model_estimation: str
        Path to the template file filled with SPM.mat information, ready to be launched
    """
    from os.path import abspath

    with open(template_file, "r") as file:
        filedata = file.read()
    # Replace by the real path to spm.mat
    filedata = filedata.replace("@SPMMAT", "'" + mat_file + "'")
    current_model_estimation = abspath("./current_model_creation.m")
    with open(current_model_estimation, "w+") as file:
        file.write(filedata)

    return current_model_estimation


def clean_spm_result_file(
    mat_file: str, template_file: str, method: str, threshold: float
) -> str:
    """Make a copy of the template file (for results) and replace @SPMMAT by the real path to SPM.mat (mat_file)

    Parameters
    ----------
    mat_file: str
        Path to the SPM.mat file of the SPM analysis
    template_file: str
        Path to the template file for getting the results of the model
    method: str
        method. In our case, "none"
    threshold: float
        cluster threshold

    Returns
    -------
    current_model_estimation: str
        Path to the template file filled with SPM.mat information, ready to be launched
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


def clean_spm_contrast_file(
    mat_file: str, template_file: str, covariates: list, class_names: list
):
    """Make a copy of the template file (for results) and replace @SPMMAT, @COVARNUMBER, @GROUP1, @GROUP2
    by the corresponding variables

    Parameters
    ----------
    mat_file: str
        Path to the SPM.mat file of the SPM analysis
    template_file: str
        Path to the template file for getting the results of the model
    covariates: list of str
        List of covariates
    class_names: list of str
        Corresponds to the 2 classes for the group comparison

    Returns
    -------
    current_model_estimation: str
        Path to the template file filled with variables, ready to be launched
    """
    from os.path import abspath

    number_of_covariates = len(covariates)

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


def copy_and_rename_spm_output_files(
    spm_mat: str,
    class_names: list,
    covariates: list,
    group_label: str,
    fwhm: int,
    measure: str,
):
    """Once analysis is done, grab all the different filenames and rename them in current directory according to class
        names

    Parameters
    ----------
    spm_mat: str
        Path to the SPM.mat file of the SPM analysis
    class_names: list of str
        Corresponds to the 2 classes for the group comparison
    covariates: list of str
        List of covariates
    group_label: str
        Name of the group label
    fwhm: int
        Fwhm in mm used
    measure: str
        Measure used

    Returns
    -------
    spmT_0001: str
        Path to t maps for the first group comparison
    spmT_0002: str
        Path to t maps for the second group comparison
    spm_figures: list
        Path to figure files
    variance_of_error: str
        Path to variance of error
    resels_per_voxels: str
        Path to resels per voxel
    mask: str
        Path to mask of included voxels
    regression_coeff: str list
        Path to regression coefficients
    contrasts: str list
        Path to weighted parameter estimation for the 2 contrasts
    """
    from os import listdir
    from os.path import abspath, dirname, isdir, isfile, join
    from shutil import copyfile

    spm_dir = dirname(spm_mat)
    if not isfile(spm_mat):
        if not isdir(spm_dir):
            raise RuntimeError("[Error] output folder " + spm_dir + " does not exist.")
        else:
            raise RuntimeError("[Error] SPM matrix " + spm_mat + " does not exist.")

    list_files = [f for f in listdir(spm_dir) if not f.startswith(".")]

    spm_figures = rename_spm_figures(spm_dir, list_files, group_label)

    spmT_0001, spmT_0002 = rename_spm_t_maps(
        spm_dir, list_files, fwhm, group_label, class_names, measure
    )

    variance_of_error = abspath(f"./group-{group_label}_VarianceError.nii")
    copyfile(abspath(join(spm_dir, "ResMS.nii")), variance_of_error)

    resels_per_voxels = abspath("./resels_per_voxel.nii")
    copyfile(abspath(join(spm_dir, "RPV.nii")), resels_per_voxels)

    mask = abspath("./included_voxel_mask.nii")
    copyfile(abspath(join(spm_dir, "mask.nii")), mask)

    regression_coeff = rename_beta_files(spm_dir, list_files, covariates, class_names)

    contrasts = rename_spm_contrast_files(
        spm_dir, list_files, group_label, class_names, measure
    )

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


def rename_spm_figures(spm_dir: str, list_files: list, group_label: str) -> list:
    """Copies and renames the figures generated by the pipeline.

    Parameters
    ----------
    spm_mat: str
        Path to the SPM.mat file of the SPM analysis
    list_files: list of str
        list of files generated by the pipeline.
    group_label: str
        Name of the group label

    Returns
    -------
    spm_figures: list
        Path to figure files
    """
    from os.path import abspath, dirname, join
    from shutil import copyfile

    figures = [abspath(join(spm_dir, f)) for f in list_files if f.endswith("png")]
    if len(figures) < 2:
        raise RuntimeError("[Error] Figures were not generated")
    fig_number = [int(f[-7:-4]) for f in figures]
    spm_figures = [
        abspath("./" + "group-" + group_label + "_report-" + str(i) + ".png")
        for i in fig_number
    ]
    for old_name, new_name in zip(figures, spm_figures):
        copyfile(old_name, new_name)
    return spm_figures


def rename_spm_t_maps(
    spm_dir: str,
    list_files: list,
    fwhm: int,
    group_label: str,
    class_names: list,
    measure: str,
) -> ty.Tuple[str]:
    """Copies and renames the spm t maps.

    Parameters
    ----------
    spm_mat: str
        Path to the SPM.mat file of the SPM analysis
    list_files: list of str
        list of files generated by the pipeline.
    fwhm: int
        Fwhm in mm used
    group_label: str
        Name of the group label
    class_names: list of str
        Corresponds to the 2 classes for the group comparison
    measure: str
        Measure used

    Returns
    -------
    spmT_0001: str
        Path to t maps for the first group comparison
    spmT_0002: str
        Path to t maps for the second group comparison
    """
    from os.path import abspath, dirname, join
    from shutil import copyfile

    spm_T = [abspath(join(spm_dir, f)) for f in list_files if f.startswith("spmT")]
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
    copyfile(join(spm_dir, "spmT_0001.nii"), spmT_0001)
    copyfile(join(spm_dir, "spmT_0002.nii"), spmT_0002)

    return spmT_0001, spmT_0002


def rename_spm_contrast_files(
    spm_dir: str, list_files: list, group_label: str, class_names: list, measure: str
) -> list:
    """Handles the contrast files.

    Parameters
    ----------
    spm_mat: str
        Path to the SPM.mat file of the SPM analysis
    list_files: list of str
        list of files generated by the pipeline.
    group_label: str
        Name of the group label
    class_names: list of str
        Corresponds to the 2 classes for the group comparison
    measure: str
        Measure used

    Returns
    -------
    contrast: list of str
        contrast
    """
    from os.path import abspath, dirname, join
    from shutil import copyfile

    con_files = [abspath(join(spm_dir, f)) for f in list_files if f.startswith("con_")]
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
    return contrasts


def rename_beta_files(
    spm_dir: str, list_files: list, covariates: list, class_names: list
) -> list:
    """Handles the beta files.

    Parameters
    ----------
    spm_mat: str
        Path to the SPM.mat file of the SPM analysis
    list_files: list of str
        list of files generated by the pipeline.
    covariates: list of str
        List of covariates
    class_names: list of str
        Corresponds to the 2 classes for the group comparison

    Returns
    -------
    regression_coeff: str list
        Path to regression coefficients
    """
    from os.path import abspath, dirname, join
    from shutil import copyfile

    betas = [abspath(join(spm_dir, f)) for f in list_files if f.startswith("beta_")]
    betas = sorted(betas)
    if len(betas) != 2 + len(covariates):
        raise RuntimeError("[Error] Not enough betas files found in output directory")
    regression_coeff_covar = [abspath(f"./{covar}.nii") for covar in covariates]
    regression_coeff = [abspath(f"./{c}.nii") for c in class_names]
    regression_coeff.extend(regression_coeff_covar)
    # Order is respected:
    for beta, reg_coeff in zip(betas, regression_coeff):
        copyfile(beta, reg_coeff)
    return regression_coeff
