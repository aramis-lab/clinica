

def clinica_file_reader(subjects,
                        sessions,
                        input_directory,
                        pattern,
                        recursive_search_max=10):
    """
    :param subjects: list of subjects
    :param sessions: list of sessions (must be same size as subjects, and must correspond )
    :param input_directory: location of the bids (or caps ?) directory
    :param pattern: define the pattern of the final file
    :param recursive_search_max: number of folder deep the function can search for the file matching the pattern
    :return: list of files respecting the subject/session order provided in input, and an error string that can have the
            following values : None (no error found) or a string describing the problem.

        raise:
            - RuntimeError if for some couple subject/session, file is not found or multiple files are found

        Examples: (path are shortened for readability)
            - You have the full name of a file:
                File orig_nu.mgz from FreeSurfer of subject sub-ADNI011S4105 session ses-M00 located in mri folder of
                FreeSurfer output :
                    clinica_file_reader(['sub-ADNI011S4105'], ['ses-M00'], caps_directory, 'orig_nu.mgz')
                    gives: ['/caps/subjects/sub-ADNI011S4105/ses-M00/t1/freesurfer_cross_sectional/sub-ADNI011S4105_ses-M00/mri/orig_nu.mgz']

            - You have a partial name of the file:
                File sub-ADNI011S4105_ses-M00_task-rest_acq-FDG_pet.nii.gz in BIDS directory. Here, filename depends on
                subject and session name :
                     clinica_file_reader(['sub-ADNI011S4105'], ['ses-M00'], bids_directory, '*fdg_pet.nii*')
                     gives: ['/bids/sub-ADNI011S4105/ses-M00/pet/sub-ADNI011S4105_ses-M00_task-rest_acq-FDG_pet.nii.gz']

            - Tricky example:
                Get the file rh.white from FreeSurfer:
                If you try:
                    clinica_file_reader(['sub-ADNI011S4105'], ['ses-M00'], caps, 'rh.white')
                        the following error will arise:
                        * More than 1 file found for pattern 'rh.white' for subject sub-Adni011S4105 and session ses-M00:
                            /caps/subjects/sub-ADNI011S4105/ses-M00/t1/freesurfer_cross_sectional/fsaverage/surf/rh.white
                            /caps/subjects/sub-ADNI011S4105/ses-M00/t1/freesurfer_cross_sectional/rh.EC_average/surf/rh.white
                            /caps/subjects/sub-ADNI011S4105/ses-M00/t1/freesurfer_cross_sectional/sub-ADNI011S4105_ses-M00/surf/rh.white
                Correct usage (wanted for pet-surface):
                    clinica_file_reader(['sub-ADNI011S4105'], ['ses-M00'], caps, 'sub-*_ses-*/surf/rh.white')
                    gives: ['/caps/subjects/sub-ADNI011S4105/ses-M00/t1/freesurfer_cross_sectional/sub-ADNI011S4105_ses-M00/surf/rh.white']

        Note:
            This function is case insensitive, meaning that the pattern argument can, for example, contain maj letter
            that do not exists in the existing file path.

    """

    from glob import glob
    from os.path import join, isdir
    from os import listdir
    from clinica.utils.io import check_bids_folder, check_caps_folder

    def insensitive_glob(pattern_glob):
        def either(c):
            return '[%s%s]' % (c.lower(), c.upper()) if c.isalpha() else c
        return glob(''.join(map(either, pattern_glob)))

    def determine_caps_or_bids(input_dir):
        """
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
                raise RuntimeError('Could not determine if ' + input_dir + ' is a CAPS or BIDS directory')

    is_bids = determine_caps_or_bids(input_directory)
    if is_bids:
        check_bids_folder(input_directory)
    else:
        check_caps_folder(input_directory)

    assert pattern[0] != '/', 'pattern argument cannot start with char : / (does not work in os.path.join function). ' \
                              + 'If you want to indicate the exact name of the file, use the format' \
                              + ' directory_name/filename.extension in the pattern argument'
    assert recursive_search_max >= 1, 'recursive_search_max argument must be >= 1'
    assert len(subjects) == len(sessions), 'Subjects and sessions must have the same length'
    rez = []
    error_encountered = []
    for sub, ses in zip(subjects, sessions):
        if is_bids:
            origin_pattern = join(input_directory, sub, ses)
        else:
            origin_pattern = join(input_directory, 'subjects', sub, ses)

        for level in range(recursive_search_max):
            current_pattern = join(origin_pattern, '*/' * level, pattern)
            current_glob_found = insensitive_glob(current_pattern)
            print('Trying ' + str(current_pattern))
            if len(current_glob_found) > 0:
                print('Found level ' + str(level))
                break

        if len(current_glob_found) > 1:
            error_str = '* More than 1 file found for pattern \'' + pattern \
                        + '\' for subject ' + sub + ' and session ' + ses + ' :\n'
            for found_file in current_glob_found:
                error_str += '\t' + found_file + '\n'
            error_encountered.append(error_str)
        elif len(current_glob_found) == 0:
            error_encountered.append('* No file found for pattern \'' + join(origin_pattern, '*/' + pattern)
                                     + '\' for subject ' + sub + ' and session ' + ses + ' until level ' + str(recursive_search_max))
        else:
            rez.append(current_glob_found[0])

    if len(error_encountered) > 0:
        error_message = 'Clinica encountered ' + str(len(error_encountered)) \
                        + ' problem(s) while getting file(s) with pattern ' + pattern + ' :\n'
        for msg in error_encountered:
                error_message += msg
    else:
        error_message = None
    return rez, error_message


if __name__ == '__main__':
    subjs = ['sub-Adni011S4105']
    sesss = ['ses-M00']
    bids = '/Users/arnaud.marcoux/CI/new_data/PETSurface/in/bids'
    file = clinica_file_reader(subjs, sesss, bids, '*fdg_pet.nii*')
    print(file)
    caps = '/Users/arnaud.marcoux/CI/new_data/PETSurface/in/caps'
    print(clinica_file_reader(subjs, sesss, caps, 'orig_nu.mgz'))
    print(clinica_file_reader(subjs, sesss, bids, '*_pet.json'))
    print(clinica_file_reader(subjs, sesss, caps, 'sub-*_ses-*/surf/rh.white'))
    print(clinica_file_reader(subjs, sesss, caps, 'sub-*_ses-*/*/lh.aparc.a2009s.annot'))
    print(clinica_file_reader(subjs, sesss, caps, 'sub-*_ses-*/*/lh.aparc.annot'))
