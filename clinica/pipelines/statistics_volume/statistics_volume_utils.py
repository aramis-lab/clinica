"""Utility functions for the pipeline statistics volume.

Public functions are called within Nipype nodes in the pipeline.
These functions must be "self-contained". In other words, they should
have type hints composed of Python's builtins only, and they should
explicitly import everything they use (i.e. other functions).

Private functions are called by public functions of this module.
They can use more advanced abstractions.
"""

import functools
import operator
import typing as ty
from pathlib import Path


def get_group_1_and_2(tsv: str, contrast: str) -> tuple:
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

    class_names = sorted(list(set(tsv[contrast])))

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
    from os import remove
    from os.path import abspath, isfile

    import pandas as pds

    from clinica.pipelines.statistics_volume.statistics_volume_utils import (  # noqa
        _convert_to_numeric,
        _create_spm_output_folder,
        _set_output_and_groups,
        _write_covariates,
    )
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

    output_folder = _create_spm_output_folder(current_model)

    # Replace string in matlab file to set the output directory, and all the scans used for group 1 and 2
    _set_output_and_groups(
        output_folder,
        current_model,
        file_list,
        idx_group1,
        idx_group2,
        filedata,
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
            current_covar_data = _convert_to_numeric(current_covar_data)
        _write_covariates(
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


def _create_spm_output_folder(current_model: str) -> str:
    """Creates the spm output folder.

    The folder is named '2_sample_t_test' and is located in the
    same folder as the provided model file.
    If the folder already exists, it is deleted and created again.

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


def _set_output_and_groups(
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
        _unravel_list_for_matlab,
    )

    filedata = filedata.replace("@OUTPUTDIR", "'" + output_folder + "'")
    filedata = filedata.replace(
        "@SCANS1",
        _unravel_list_for_matlab([file_list[i] for i in idx_group1]),
    )
    filedata = filedata.replace(
        "@SCANS2",
        _unravel_list_for_matlab([file_list[i] for i in idx_group2]),
    )
    with open(current_model, "w+") as file:
        file.write(filedata)


def _unravel_list_for_matlab(my_list: ty.List[str]) -> str:
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


def _convert_to_numeric(data: ty.List[str]) -> ty.List[float]:
    """Convert the string values from the input list to numeric values.

    If the values can be casted to floats, then outputs a list of floats.
    If the values are categorical, then performs encoding of the categories.

    Parameters
    ----------
    data: list of str
        List of strings to be converted.

    Returns
    -------
    list :
        List of floats resulting from conversion.
    """
    import numpy as np

    from clinica.pipelines.statistics_volume.statistics_volume_utils import (  # noqa
        _is_number,
    )

    temp_data = [elem.replace(",", ".") for elem in data]
    if all(_is_number(elem) for elem in temp_data):
        return [float(elem) for elem in temp_data]
    # categorical variables (like Male; Female; M, F etc...)
    unique_values = list(np.unique(np.array(data)))
    return [unique_values.index(elem) for elem in data]


def _is_number(s: str) -> bool:
    """Returns True is the input can be converted to float, False otherwise."""
    try:
        # Try the conversion, if it is not possible, error will be raised
        float(s)
        return True
    except ValueError:
        return False
    except TypeError:
        return False


def _write_covariates(
    covariates: list,
    idx_group1: list,
    idx_group2: list,
    model: str,
    covar: int,
    covar_number: int,
) -> None:
    """Concatenate covariates for group1 and group2 and write them to disk.

    Parameters
    ----------
    covariates: list
        List of covariates to be concatenated and written.
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
    from clinica.pipelines.statistics_volume.statistics_volume_utils import (  # noqa
        _write_covariate_lines,
    )

    current_covar_data_group1 = [covariates[i] for i in idx_group1]
    current_covar_data_group2 = [covariates[i] for i in idx_group2]
    concatenated_covariates = current_covar_data_group1 + current_covar_data_group2
    _write_covariate_lines(model, covar_number, covar, concatenated_covariates)


def _write_covariate_lines(
    m_file_to_write_in: str,
    covar_number: int,
    covar_name: str,
    covar_values: list,
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
    from pathlib import Path

    from clinica.pipelines.statistics_volume.statistics_volume_utils import (  # noqa
        _delete_last_line,
        _run_matlab_script_with_matlab,
        _run_matlab_script_with_spm_standalone,
    )
    from clinica.utils.spm import spm_standalone_is_available

    m_file = Path(m_file)
    if not m_file.exists():
        raise FileNotFoundError(f"[Error] File {m_file} does not exist")
    if m_file.suffix != ".m":
        raise ValueError(
            f"[Error] {m_file} is not a Matlab file (extension must be .m)"
        )
    if spm_standalone_is_available():
        _delete_last_line(m_file)
        _run_matlab_script_with_spm_standalone(m_file)
    else:
        _run_matlab_script_with_matlab(str(m_file))
    output_mat_file = (m_file.parent.parent / "2_sample_t_test" / "SPM.mat").resolve()
    if not output_mat_file.is_file():
        raise RuntimeError(f"Output matrix {output_mat_file} was not produced")
    return str(output_mat_file)


def _delete_last_line(filename: Path) -> None:
    """Removes the last line of the provided file.

    Used to remove the call to spm jobman if the MATLAB file
    is used with SPM standalone.

    Parameters
    ----------
    filename: Path
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


def _run_matlab_script_with_spm_standalone(m_file: Path) -> None:
    """Runs a matlab script spm_standalone

    Parameters
    ----------
    m_file: str
        Path to the script
    """
    from os import system

    system(_get_matlab_standalone_command(m_file))


def _get_matlab_standalone_command(m_file: Path) -> str:
    """Get the matlab standalone command to be run depending on the user's OS.

    Notes
    -----
    On MAC systems, SPM standalone must be run directly from its root folder.

    Raises
    ------
    SystemError
        If this is run on unsupported platforms.
    """
    if _is_running_on_mac():
        return f"cd $SPMSTANDALONE_HOME && ./run_spm12.sh $MCR_HOME batch {m_file}"
    if _is_running_on_linux():
        return f"$SPMSTANDALONE_HOME/run_spm12.sh $MCR_HOME batch {m_file}"
    raise SystemError("Clinica only support Mac OS and Linux")


def _is_running_on_os(os: str) -> bool:
    import platform

    return platform.system().lower().startswith(os)


_is_running_on_mac = functools.partial(_is_running_on_os, os="darwin")
_is_running_on_linux = functools.partial(_is_running_on_os, os="linux")


def _run_matlab_script_with_matlab(m_file: str) -> None:
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
    matlab.inputs.paths = Path(m_file).parent
    matlab.inputs.script = Path(m_file).stem
    matlab.inputs.single_comp_thread = False
    matlab.inputs.logfile = abspath("./matlab_output.log")
    matlab.run()


def clean_template_file(mat_file: str, template_file: str) -> str:
    """Make a copy of the template file (for estimation) and replace
    @SPMMAT by the real path to SPM.mat (mat_file).

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
    mat_file: str,
    template_file: str,
    method: str,
    threshold: float,
) -> str:
    """Make a copy of the template file (for results) and replace
    @SPMMAT by the real path to SPM.mat (mat_file).

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
    mat_file: str,
    template_file: str,
    covariates: list,
    class_names: list,
):
    """Make a copy of the template file (for results) and replace
    @SPMMAT, @COVARNUMBER, @GROUP1, @GROUP2 by the corresponding variables.

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
    output_dir: str = ".",
) -> list:
    """Once analysis is done, copy and rename (according to class names)
    the different filenames in the provided output directory.

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
    output_dir : str, optional
        Output folder where the copied figures should be written.
        Default to current folder.

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
    from clinica.pipelines.statistics_volume.statistics_volume_utils import (  # noqa
        _rename_beta_files,
        _rename_spm_contrast_files,
        _rename_spm_figures,
        _rename_spm_t_maps,
    )

    spm_dir, output_dir = _check_spm_and_output_dir(spm_mat, output_dir)
    parameters = {
        "group_label": group_label,
        "measure": measure,
        "class_names": class_names,
        "fwhm": fwhm,
        "covariates": covariates,
    }
    rename_functions = (
        _rename_spm_t_maps,
        _rename_spm_figures,
        _rename_other_spm_files,
        _rename_beta_files,
        _rename_spm_contrast_files,
    )
    return [
        func(
            spm_dir=spm_dir,
            output_dir=output_dir,
            parameters_for_new_filename_construction=parameters,
        )
        for func in rename_functions
    ]


def _rename_spm_t_maps(
    spm_dir: Path,
    output_dir: Path,
    parameters_for_new_filename_construction: dict,
) -> ty.List[str]:
    from clinica.pipelines.statistics_volume.statistics_volume_utils import (  # noqa
        _build_t_map_filenames,
        _get_spm_t_maps,
    )

    return _rename_spm_files(
        spm_dir,
        output_dir,
        spm_filename_getter=_get_spm_t_maps,
        new_filename_constructor=_build_t_map_filenames,
        parameters_for_new_filename_construction=parameters_for_new_filename_construction,
    )


def _rename_spm_files(
    spm_dir: Path,
    output_dir: Path,
    spm_filename_getter,
    new_filename_constructor,
    parameters_for_new_filename_construction: dict,
) -> ty.List[str]:
    """bla"""
    kwargs = dict()
    if spm_filename_getter.__name__ == "_get_spm_beta_files":
        kwargs["expected_number_of_files"] = 2 + len(
            parameters_for_new_filename_construction["covariates"]
        )
    existing_spm_filenames = spm_filename_getter(spm_dir, **kwargs)
    new_filenames = new_filename_constructor(
        **parameters_for_new_filename_construction,
        existing_spm_filenames=existing_spm_filenames,
    )
    spm_to_clinica_filename_mapping = {
        k: v for k, v in zip(existing_spm_filenames, new_filenames)
    }
    _copy_spm_to_clinica_files(spm_dir, output_dir, spm_to_clinica_filename_mapping)

    return [str(output_dir / f) for f in new_filenames]


def _copy_spm_to_clinica_files(
    spm_dir: Path, output_dir: Path, name_mapping: ty.Dict[str, str]
) -> None:
    _copy_files(
        {str(spm_dir / k): str(output_dir / v) for k, v in name_mapping.items()}
    )


def _copy_files(source_to_destination_mapping: ty.Dict[str, str]) -> None:
    """Copy the source files in the mapping keys to the destinations in the values."""
    from shutil import copyfile

    for source, destination in source_to_destination_mapping.items():
        copyfile(source, destination)


def _rename_spm_contrast_files(
    spm_dir: Path,
    output_dir: Path,
    parameters_for_new_filename_construction: dict,
) -> ty.List[str]:
    from clinica.pipelines.statistics_volume.statistics_volume_utils import (  # noqa
        _get_new_spm_contrast_files,
        _get_spm_contrast_files,
    )

    return _rename_spm_files(
        spm_dir,
        output_dir,
        spm_filename_getter=_get_spm_contrast_files,
        new_filename_constructor=_get_new_spm_contrast_files,
        parameters_for_new_filename_construction=parameters_for_new_filename_construction,
    )


def _rename_beta_files(
    spm_dir: Path,
    output_dir: Path,
    parameters_for_new_filename_construction: dict,
) -> ty.List[str]:
    from clinica.pipelines.statistics_volume.statistics_volume_utils import (  # noqa
        _get_new_spm_beta_files,
        _get_spm_beta_files,
    )

    return _rename_spm_files(
        spm_dir,
        output_dir,
        spm_filename_getter=_get_spm_beta_files,
        new_filename_constructor=_get_new_spm_beta_files,
        parameters_for_new_filename_construction=parameters_for_new_filename_construction,
    )


def _rename_spm_figures(
    spm_dir: Path,
    output_dir: Path,
    parameters_for_new_filename_construction: dict,
) -> ty.List[str]:
    from clinica.pipelines.statistics_volume.statistics_volume_utils import (  # noqa
        _get_new_spm_figures,
        _get_spm_figures,
    )

    return _rename_spm_files(
        spm_dir,
        output_dir,
        spm_filename_getter=_get_spm_figures,
        new_filename_constructor=_get_new_spm_figures,
        parameters_for_new_filename_construction=parameters_for_new_filename_construction,
    )


def _rename_other_spm_files(
    spm_dir: Path,
    output_dir: Path,
    parameters_for_new_filename_construction: dict,
) -> ty.List[str]:
    from clinica.pipelines.statistics_volume.statistics_volume_utils import (  # noqa
        _get_new_other_spm_files,
        _get_other_spm_files,
    )

    return _rename_spm_files(
        spm_dir,
        output_dir,
        spm_filename_getter=_get_other_spm_files,
        new_filename_constructor=_get_new_other_spm_files,
        parameters_for_new_filename_construction=parameters_for_new_filename_construction,
    )


def _get_spm_t_maps(spm_dir: Path) -> ty.List[str]:
    return _get_spm_files(spm_dir, pattern="spmT", file_type="t-map(s)")


def _get_other_spm_files(spm_dir: Path, **kwargs) -> ty.List[str]:
    return ["ResMS.nii", "RPV.nii", "mask.nii"]


def _get_spm_figures(spm_dir: Path) -> ty.List[str]:
    import operator

    return _get_spm_files(
        spm_dir,
        pattern="png",
        grab_by_prefix=False,
        file_type="figures",
        op=operator.ge,
    )


def _get_spm_beta_files(spm_dir: Path, expected_number_of_files: int) -> ty.List[str]:
    return _get_spm_files(
        spm_dir,
        pattern="beta_",
        expected_number_of_files=expected_number_of_files,
        file_type="betas",
    )


def _get_spm_contrast_files(spm_dir: Path) -> ty.List[str]:
    return _get_spm_files(spm_dir, pattern="con_", file_type="contrast")


def _get_spm_files(
    spm_dir: Path,
    pattern: str,
    file_type: str,
    grab_by_prefix: bool = True,
    expected_number_of_files: int = 2,
    op=operator.eq,
) -> ty.List[str]:
    """Add docstring here."""
    spm_files = [
        f.name
        for f in spm_dir.iterdir()
        if _is_valid_filename(f.name, pattern, grab_by_prefix)
    ]
    if not op(len(spm_files), expected_number_of_files):
        raise RuntimeError(
            f"[Error] Wrong number of SPM {file_type} files. "
            f"Expected {expected_number_of_files}, got {len(spm_files)}."
        )

    return spm_files


def _is_valid_filename(filename: str, pattern: str, grab_by_prefix: bool) -> bool:
    """Filename is valid if it doesn't start with '.' and if it has the provided
    pattern as a prefix or as a suffix.
    """
    pattern_condition = (
        filename.startswith(pattern) if grab_by_prefix else filename.endswith(pattern)
    )

    return not filename.startswith(".") and pattern_condition


def _get_new_other_spm_files(group_label: str, **kwargs):
    return [
        f"group-{group_label}_VarianceError.nii",
        "resels_per_voxel.nii",
        "included_voxel_mask.nii",
    ]


def _get_new_spm_figures(
    group_label: str, existing_spm_filenames: ty.List[str], **kwargs
) -> ty.List[str]:
    fig_number = [
        int(f[-7:-4]) for f in existing_spm_filenames
    ]  # assumes number is encoded on 3 digits
    return [f"group-{group_label}_report-{i}.png" for i in fig_number]


def _check_spm_and_output_dir(
    spm_mat: str, output_dir: ty.Optional[str] = None
) -> ty.Tuple[Path, Path]:
    from pathlib import Path

    spm_mat = Path(spm_mat)
    spm_dir = spm_mat.parent

    if not spm_mat.is_file():
        if not spm_dir.is_dir():
            raise RuntimeError(f"[Error] output folder {spm_dir} does not exist.")
        raise RuntimeError(f"[Error] SPM matrix {spm_mat} does not exist.")

    output_dir = Path(output_dir) if output_dir else Path(".")

    return spm_dir.resolve(), output_dir.resolve()


def _build_t_map_filenames(
    class_names: ty.List[str],
    group_label: str,
    measure: str,
    fwhm: int,
    **kwargs,
) -> ty.List[str]:
    """Build the T-map filenames from the class names, group label, measure, and fwhm."""
    fwhm_string = f"_fwhm-{int(fwhm)}" if fwhm else ""

    return [
        f"group-{group_label}_{contrast}_measure-{measure}{fwhm_string}_TStatistics.nii"
        for contrast in _build_contrasts_from_class_names(class_names)
    ]


def _build_contrasts_from_class_names(class_names: ty.List[str]) -> ty.List[str]:
    """Build the contrast strings from the class names."""
    return [
        f"{class_names[0]}-lt-{class_names[1]}",
        f"{class_names[1]}-lt-{class_names[0]}",
    ]


def _get_new_spm_beta_files(
    class_names: ty.List[str], covariates: ty.List[str], **kwargs
) -> ty.List[str]:
    return [f"{c}.nii" for c in class_names] + [f"{covar}.nii" for covar in covariates]


def _get_new_spm_contrast_files(
    group_label: str, class_names: ty.List[str], measure: str, **kwargs
) -> ty.List[str]:
    return [
        f"group-{group_label}_{contrast}_measure-{measure}_contrast.nii"
        for contrast in _build_contrasts_from_class_names(class_names)
    ]
