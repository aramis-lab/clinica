# coding: utf8


def get_group_1_and_2(csv, contrast):
    """
    @TODO write doc
    Args:
        csv_file:
        contrast:

    Returns:

    """
    import pandas as pds
    from clinica.utils.exceptions import ClinicaException

    csv = pds.read_csv(csv, sep='\t')
    columns = list(csv.columns)

    if contrast not in columns:
        raise ClinicaException(contrast + ' is not present in ' + csv)

    values_of_contrast = list(set(csv[contrast]))
    if len(values_of_contrast) != 2:
        raise ClinicaException('It must exist 2 classes for the contrast category. Here Clinica found: '
                               + str(values_of_contrast))
    first_group_idx = [i for i, label in enumerate(list(csv[contrast])) if label == values_of_contrast[0]]
    second_group_idx = [i for i, label in enumerate(list(csv[contrast])) if label == values_of_contrast[1]]

    return first_group_idx, second_group_idx


def model_creation(csv, contrast, idx_group1, idx_group2, file_list, template_file):
    """

    Args:
        csv:
        contrast:
        idx_group1:
        idx_group2:
        file_list:
        template_file:

    Returns:

    """
    from os.path import join, dirname, isfile, abspath
    from clinica.utils.exceptions import ClinicaException
    from numbers import Number
    from os import remove, mkdir
    import numpy as np
    import pandas as pds
    import clinica.pipelines.statistics_volume.statistics_volume_utils as utls

    # Get template for model creation
    if not isfile(template_file):
        raise RuntimeError('[Error] ' + template_file + ' could not be found !')

    # Create current model filename
    current_model = abspath('./current_model_creation.m')

    if isfile(current_model):
        remove(current_model)

    # Read template
    with open(template_file, 'r') as file:
        filedata = file.read()
    output_folder = abspath(join(dirname(current_model), '..', '2_sample_t_test'))
    mkdir(output_folder)
    # Replace our string
    filedata = filedata.replace('@OUTPUTDIR', '\'' + output_folder + '\'')
    filedata = filedata.replace('@SCANS1', utls.unravel_list_for_matlab([f for i, f in enumerate(file_list) if i in idx_group1]))
    filedata = filedata.replace('@SCANS2', utls.unravel_list_for_matlab([f for i, f in enumerate(file_list) if i in idx_group2]))

    with open(current_model, 'w+') as file:
        file.write(filedata)

    # Add our covariates
    csv = pds.read_csv(csv, sep='\t')
    columns_stripped = [elem.strip(' ') for elem in list(csv.columns)]
    if columns_stripped != list(csv.columns):
        raise ClinicaException('[Error] Check the column of your tsv file ' + csv
                               + 'Whitespace in the column names can cause errors')
    covariables = [elem for elem in columns_stripped if elem not in ['participant_id', 'session_id', contrast]]
    for covar_number, covar in enumerate(covariables, start=1):
        current_covar_data = list(csv[covar])
        if isinstance(current_covar_data[0], str):
            # Transform data
            temp_data = [elem.replace(',', '.') for elem in current_covar_data]
            if all(utls.is_number(elem) for elem in temp_data):
                current_covar_data = [float(elem) for elem in temp_data]
            else:
                # categorical variables (like Male; Female; M, F etc...)
                unique_values = list(np.unique(np.array(current_covar_data)))
                current_covar_data = [unique_values.index(elem) for elem in current_covar_data]

        elif isinstance(current_covar_data[0], Number):
            # Do nothing
            pass

        current_covar_data_group1 = [elem for i, elem in enumerate(current_covar_data) if i in idx_group1]
        current_covar_data_group2 = [elem for i, elem in enumerate(current_covar_data) if i in idx_group2]
        covar_data_concatenated = current_covar_data_group1 + current_covar_data_group2
        utls.write_covariable_lines(current_model, covar_number, covar, covar_data_concatenated)

    # Tell matlab to run the script at the end
    with open(current_model, 'a') as file:
        file.write('spm_jobman(\'run\', matlabbatch)')
    return current_model
    

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def unravel_list_for_matlab(my_list):
    """

    Args:
        my_list:

    Returns:

    """
    result = '\',\''.join(my_list)
    result = '\'' + result + '\''
    return result


def write_covariable_lines(m_file_to_write_in, covar_number, covar_name, covar_values):
    """

    Args:
        m_file_to_write_in:
        covar_number: must start at 1
        covar_name:
        covar_values:

    Returns:

    """
    from os.path import isfile

    assert isfile(m_file_to_write_in), 'Could not find file ' + m_file_to_write_in
    covar_values_string = [str(elem) for elem in covar_values]
    covar_values_for_script = '[' + ' '.join(covar_values_string) + ']'

    m_file = open(m_file_to_write_in, 'a')
    m_file.write('matlabbatch{1}.spm.stats.factorial_design.cov(' + str(covar_number) + ').c = ' + covar_values_for_script + ';\n')
    m_file.write('matlabbatch{1}.spm.stats.factorial_design.cov(' + str(covar_number) + ').cname = \'' + covar_name + '\';\n')
    m_file.write('matlabbatch{1}.spm.stats.factorial_design.cov(' + str(covar_number) + ').iCFI = 1;\n')
    m_file.write('matlabbatch{1}.spm.stats.factorial_design.cov(' + str(covar_number) + ').iCC = 1;\n')
    m_file.close()


def run_m_script(m_file):
    from os.path import isfile, dirname, basename, abspath, join
    from os import system
    import clinica.pipelines.statistics_volume.statistics_volume_utils as utls
    from nipype.interfaces.matlab import MatlabCommand, get_matlab_command
    import platform
    from clinica.utils.stream import cprint

    assert isinstance(m_file, str), '[Error] Argument must be a string'
    if not isfile(m_file):
        raise FileNotFoundError('[Error] File ' + m_file + 'does not exist')
    assert m_file[-2:] == '.m', '[Error] ' + m_file + ' is not a Matlab file (extension must be .m)'

    # Generate command line to run
    if utls.use_spm_standalone():
        cprint('USING SPM STANDALONE')
        utls.delete_last_line(m_file)
        # SPM standalone must be run directly from its root folder
        if platform.system().lower().startswith('darwin'):
            # Mac OS
            cmdline = 'cd $SPMSTANDALONE_HOME && ./run_spm12.sh $MCR_HOME batch ' + m_file
        elif platform.system().lower().startswith('linux'):
            # Linux OS
            cmdline = '$SPMSTANDALONE_HOME/run_spm12.sh $MCR_HOME batch ' + m_file
        else:
            raise SystemError('Clinica only support Mac OS and Linux')
        system(cmdline)
    else:
        cprint('USING CLASSICAL SPM')
        MatlabCommand.set_default_matlab_cmd(get_matlab_command())
        matlab = MatlabCommand()
        if platform.system().lower().startswith('linux'):
            matlab.inputs.args = '-nosoftwareopengl'
        matlab.inputs.paths = dirname(m_file)
        matlab.inputs.script = basename(m_file)[:-2]
        matlab.inputs.single_comp_thread = False
        matlab.inputs.logfile = abspath('./matlab_output.log')
        matlab.run()
    output_mat_file = abspath(join(dirname(m_file), '..', '2_sample_t_test', 'SPM.mat'))
    if not isfile(output_mat_file):
        raise RuntimeError('Output matrix ' + output_mat_file + ' was not produced')
    return output_mat_file


def use_spm_standalone():
    import os
    from os.path import isdir, expandvars

    use_spm_standalone = False
    if all(elem in os.environ.keys() for elem in ['SPMSTANDALONE_HOME', 'MCR_HOME']):
        if isdir(expandvars('$SPMSTANDALONE_HOME')) and isdir(expandvars('$MCR_HOME')):
            use_spm_standalone = True
        else:
            raise FileNotFoundError('[Error] $SPMSTANDALONE_HOME and $MCR_HOME are defined, but linked to non existent folder')
    return use_spm_standalone


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

    Args:
        mat_file:
        template_file:

    Returns:

    """
    from os.path import abspath

    # Read template
    with open(template_file, 'r') as file:
        filedata = file.read()
    # Replace by the real path to spm.mat
    filedata = filedata.replace('@SPMMAT', '\'' + mat_file + '\'')
    current_model_estimation = abspath('./current_model_creation.m')
    with open(current_model_estimation, 'w+') as file:
        file.write(filedata)

    return current_model_estimation


def results(mat_file, template_file):
    """

    Args:
        mat_file:
        template_file:

    Returns:

    """
    from os.path import abspath

    # Read template
    with open(template_file, 'r') as file:
        filedata = file.read()
    # Replace by the real path to spm.mat
    filedata = filedata.replace('@SPMMAT', '\'' + mat_file + '\'')
    current_model_result = abspath('./current_model_result.m')
    with open(current_model_result, 'w+') as file:
        file.write(filedata)

    return current_model_result
