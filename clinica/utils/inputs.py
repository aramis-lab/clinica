def insensitive_glob(pattern_glob):
    """
    This function is the glob.glob() function that is insensitive to the case
    :param pattern_glob: sensitive-to-the-case pattern
    :return: insensitive-to-the-case pattern
    """
    from glob import glob

    def either(c):
        return '[%s%s]' % (c.lower(), c.upper()) if c.isalpha() else c

    return glob(''.join(map(either, pattern_glob)))


def clinica_file_reader(subjects,
                        sessions,
                        input_directory,
                        information,
                        recursive_search_max=10):
    """
    This function grabs files relative to a subject and session list according to a glob pattern (using *)
    Args:
        subjects: list of subjects
        sessions: list of sessions (must be same size as subjects, and must correspond )
        input_directory: location of the bids (or caps ?) directory
        information: dictionnary containg all the relevant information to look for the files. Dict must contains the
                     following keys : pattern, description, needed_pipeline
                             pattern: define the pattern of the final file
                             description: string to describe what the file is
                             needed_pipeline (optional): string describing the pipeline needed to obtain the file beforehand
        recursive_search_max: number of folder deep the function can search for the file matching the pattern

    Returns:
         list of files respecting the subject/session order provided in input, and an error string that can have the
         following values : None (no error found) or a string describing the problem.

        raise:
            - Nothing, we prefer returning a string with the problem written in it, so that the rest of Clinica can
            handle the problem properly

        Examples: (path are shortened for readability)
            - You have the full name of a file:
                File orig_nu.mgz from FreeSurfer of subject sub-ADNI011S4105 session ses-M00 located in mri folder of
                FreeSurfer output :
                    clinica_file_reader(['sub-ADNI011S4105'],
                                        ['ses-M00'],
                                        caps_directory,
                                        {'pattern': 'orig_nu.mgz',
                                         'description': 'freesurfer file orig_nu.mgz',
                                         'needed_pipeline': 't1-freesurfer'})
                    gives: ['/caps/subjects/sub-ADNI011S4105/ses-M00/t1/freesurfer_cross_sectional/sub-ADNI011S4105_ses-M00/mri/orig_nu.mgz']

            - You have a partial name of the file:
                File sub-ADNI011S4105_ses-M00_task-rest_acq-FDG_pet.nii.gz in BIDS directory. Here, filename depends on
                subject and session name :
                     clinica_file_reader(['sub-ADNI011S4105'],
                                         ['ses-M00'],
                                         bids_directory,
                                         {'pattern': '*fdg_pet.nii*',
                                          'description': 'FDG PET data'})
                     gives: ['/bids/sub-ADNI011S4105/ses-M00/pet/sub-ADNI011S4105_ses-M00_task-rest_acq-FDG_pet.nii.gz']

            - Tricky example:
                Get the file rh.white from FreeSurfer:
                If you try:
                    clinica_file_reader(['sub-ADNI011S4105'],
                                        ['ses-M00'],
                                        caps,
                                        {'pattern': 'rh.white',
                                         'description': 'right hemisphere of outter cortical surface.',
                                         'needed_pipeline': 't1-freesurfer'})
                        the following error will arise:
                        * More than 1 file found::
                            /caps/subjects/sub-ADNI011S4105/ses-M00/t1/freesurfer_cross_sectional/fsaverage/surf/rh.white
                            /caps/subjects/sub-ADNI011S4105/ses-M00/t1/freesurfer_cross_sectional/rh.EC_average/surf/rh.white
                            /caps/subjects/sub-ADNI011S4105/ses-M00/t1/freesurfer_cross_sectional/sub-ADNI011S4105_ses-M00/surf/rh.white
                Correct usage (wanted for pet-surface): pattern string must be 'sub-*_ses-*/surf/rh.white')
                    It then gives: ['/caps/subjects/sub-ADNI011S4105/ses-M00/t1/freesurfer_cross_sectional/sub-ADNI011S4105_ses-M00/surf/rh.white']

        Note:
            This function is case insensitive, meaning that the pattern argument can, for example, contain maj letter
            that do not exists in the existing file path.

    """

    from os.path import join, isdir
    from os import listdir
    from clinica.utils.io import check_bids_folder, check_caps_folder
    from colorama import Fore

    def determine_caps_or_bids(input_dir):
        """
        Determines if the input is a CAPS or a BIDS folder
        :param input_dir: input folder
        :return: True if input_dir is a bids, False if input_dir is a CAPS

            raise :
                RuntimeError if function could not determine if BIDS or CAPS or whatever else
        """
        if isdir(join(input_directory, 'subjects')):
            if len([sub for sub in listdir(join(input_dir, 'subjects')) if (sub.startswith('sub-') and isdir(join(input_dir, 'subjects', sub)))]) > 0 or isdir(join(input_dir, 'groups')):
                return False
            else:
                raise RuntimeError('Could not determine if ' + input_dir + ' is a CAPS or BIDS directory')

        else:
            if len([sub for sub in listdir(input_dir) if (sub.startswith('sub-') and isdir(join(input_dir, sub)))]) > 0:
                return True
            else:
                if isdir(join(input_dir, 'groups')):
                    return False
                else:
                    raise RuntimeError('Could not determine if ' + input_dir + ' is a CAPS or BIDS directory')

    assert isinstance(information, dict), 'A dict must be provided for the argmuent \'dict\''
    assert all(elem in information.keys() for elem in ['pattern', 'description']), '\'information\' must contain the keys \'pattern\' and \'description'
    assert all(elem in ['pattern', 'description', 'needed_pipeline'] for elem in information.keys()), '\'information\' can only contain the keys \'pattern\', \'description\' and \'needed_pipeline\''

    pattern = information['pattern']
    is_bids = determine_caps_or_bids(input_directory)

    if is_bids:
        check_bids_folder(input_directory)
    else:
        check_caps_folder(input_directory)

    # Some check on the formatting on the data
    assert pattern[0] != '/', 'pattern argument cannot start with char : / (does not work in os.path.join function). ' \
                              + 'If you want to indicate the exact name of the file, use the format' \
                              + ' directory_name/filename.extension or filename.extensionin the pattern argument'
    assert recursive_search_max >= 1, 'recursive_search_max argument must be >= 1'
    assert len(subjects) == len(sessions), 'Subjects and sessions must have the same length'
    if len(subjects) == 0:
        raise RuntimeError('No subjects and sessions provided.')

    # rez is the list containing the results
    rez = []
    # error is the list of the errors that happen during the whole process
    error_encountered = []
    for sub, ses in zip(subjects, sessions):
        if is_bids:
            origin_pattern = join(input_directory, sub, ses)
        else:
            origin_pattern = join(input_directory, 'subjects', sub, ses)

        # This search at each level if the file is found. We stop when a result is found or when the
        # maximum level of depth is reached
        for level in range(recursive_search_max):
            current_pattern = join(origin_pattern, '*/' * level, pattern)
            current_glob_found = insensitive_glob(current_pattern)
            # print('Trying ' + str(current_pattern))
            if len(current_glob_found) > 0:
                # print('Found level ' + str(level))
                break

        # Error handling if more than 1 file are found, or when no file is found
        if len(current_glob_found) > 1:
            error_str = '\t*' + Fore.BLUE + ' (' + sub + ' | ' + ses + ') ' + Fore.RESET + ': More than 1 file found:\n'
            for found_file in current_glob_found:
                error_str += '\t\t' + found_file + '\n'
            error_encountered.append(error_str)
        elif len(current_glob_found) == 0:
            error_encountered.append('\t*' + Fore.BLUE + ' (' + sub + ' | ' + ses + ') ' + Fore.RESET + ': No file found')
        # Otherwise the file found is added to the result
        else:
            rez.append(current_glob_found[0])

    # We do not raise an error, so that the developper can gather all the problems before Clinica crashes
    if len(error_encountered) > 0:
        error_message = Fore.RED + '\n[Error] Clinica encountered ' + str(len(error_encountered)) \
                        + ' problem(s) while getting ' + information['description'] + ':\n' + Fore.RESET
        if 'needed_pipeline' in information.keys():
            if information['needed_pipeline']:
                error_message += Fore.YELLOW + 'Please note that the following clinica pipeline(s) must have run ' \
                                 'to obtain these files: ' + information['needed_pipeline'] + Fore.RESET + '\n'
        for msg in error_encountered:
                error_message += msg
    else:
        error_message = None
    return rez, error_message


def clinica_group_reader(caps_directory, information, recursive_search_max=10):
    """
    This function grabs files relative to a group, according to a glob pattern (using *)
    Args:
        caps_directory: input caps directory
        information: dictionnary containg all the relevant information to look for the files. Dict must contains the
                     following keys : pattern, description, needed_pipeline
                             pattern: define the pattern of the final file
                             description: string to describe what the file is
                             needed_pipeline (optional): string describing the pipeline needed to obtain the file beforehand
        recursive_search_max: number of folder deep the function can search for the file matching the pattern

    Returns:
          list of files and an error string that can have the following values : None (no error found) or
          a string describing the problem.
    """
    from clinica.utils.io import check_caps_folder
    from os.path import join
    from colorama import Fore

    assert isinstance(information, dict), 'A dict must be provided for the argmuent \'dict\''
    assert all(elem in information.keys()
               for elem in ['pattern', 'description', 'needed_pipeline']), '\'information\' must contain the keys \'pattern\', \'description\', \'needed_pipeline\''

    pattern = information['pattern']
    # Some check on the formatting on the data
    assert pattern[0] != '/', 'pattern argument cannot start with char : / (does not work in os.path.join function). ' \
                              + 'If you want to indicate the exact name of the file, use the format' \
                              + ' directory_name/filename.extension or filename.extensionin the pattern argument'
    assert recursive_search_max >= 1, 'recursive_search_max argument must be >= 1'

    check_caps_folder(caps_directory)

    error_string = None
    # This search at each level if the file is found. We stop when a result is found or when the
    # maximum level of depth is reached
    for level in range(recursive_search_max):
        current_pattern = join(caps_directory, '*/' * level, pattern)
        current_glob_found = insensitive_glob(current_pattern)
        # print('Trying ' + str(current_pattern))
        if len(current_glob_found) > 0:
            break
    if len(current_glob_found) == 0:
        error_string = Fore.RED + '\n[Error] Clinica encountered a problem while getting ' + information['description'] \
                       + '.' + Fore.RESET + '\n\tCAPS directory: ' + caps_directory + '\n' + Fore.YELLOW \
                       + 'Please note that the following clinica pipeline(s) must have run to obtain these files: ' \
                       + information['needed_pipeline'] + Fore.RESET + '\n'
    return current_glob_found, error_string


if __name__ == '__main__':
    # subjs = ['sub-Adni011S4105']
    # sesss = ['ses-M00']
    # bids = '/Users/arnaud.marcoux/CI/new_data/PETSurface/in/bids'
    # file, err = clinica_file_reader(subjs,
    #                                 sesss,
    #                                 bids,
    #                                 {'pattern': '*fdg_pet.nii*',
    #                                  'description': 'FDG PET file'})
    # caps = '/Users/arnaud.marcoux/CI/new_data/PETSurface/in/caps'
    # file, err = clinica_file_reader(subjs,
    #                                 sesss,
    #                                 caps,
    #                                 {'pattern': 'surf/rh.white',
    #                                  'description': 'right hemisphere of outter cortical surface.',
    #                                  'needed_pipeline': 't1-freesurfer'})
    # print(err)
    """
    print(file)
    caps = '/Users/arnaud.marcoux/CI/new_data/PETSurface/in/caps'
    print(clinica_file_reader(subjs, sesss, caps, 'orig_nu.mgz'))
    print(clinica_file_reader(subjs, sesss, bids, '*_pet.json'))
    print(clinica_file_reader(subjs, sesss, caps, 'sub-*_ses-*/surf/rh.white'))
    """
    g_id = 'UnitTest'
    a, st = clinica_group_reader('/Users/arnaud.marcoux/CI/new_data/T1VolumeExistingTemplate/in/caps',
                                 {'pattern': 'group-' + g_id + '*_iteration-*_templaAte.nii*',
                                  'description': 'template file for ' + g_id,
                                  'needed_pipeline': 't1-volume-create_dartel or t1-volume'})
    # a, st = clinica_group_reader('/Users/arnaud.marcoux/CI/new_data/T1VolumeExistingTemplate/in/caps',
    #                             'group-' + g_id + '_template.nii*',
    #                             recursive_search_max=10)
    print(st)
