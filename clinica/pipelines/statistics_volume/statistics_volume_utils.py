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

    Returns:

    """
    from os.path import join, dirname, isfile, abspath
    from clinica.utils.exceptions import ClinicaException
    from numbers import Number
    from os import remove
    import numpy as np
    import pandas as pds
    import clinica.pipelines.statistics_volume.statistics_volume_utils as utls
    from os import system

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

    # Replace our string
    filedata = filedata.replace('@OUTPUTDIR', '\'' + dirname(current_model) + '\'')
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
