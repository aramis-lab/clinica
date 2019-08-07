# coding: utf8


def zip_nii(in_file, same_dir=False):
    from os import getcwd
    from os.path import abspath, join
    import gzip
    import shutil
    from nipype.utils.filemanip import split_filename
    from traits.trait_base import _Undefined

    if (in_file is None) or isinstance(in_file, _Undefined):
        return None

    if not isinstance(in_file, str):  # type(in_file) is list:
        return [zip_nii(f, same_dir) for f in in_file]

    orig_dir, base, ext = split_filename(str(in_file))

    # Already compressed
    if ext[-3:].lower() == ".gz":
        return in_file
    # Not compressed

    if same_dir:
        out_file = abspath(join(orig_dir, base + ext + '.gz'))
    else:
        out_file = abspath(join(getcwd(), base + ext + '.gz'))

    with open(in_file, 'rb') as f_in, gzip.open(out_file, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)

    return out_file


def unzip_nii(in_file):
    from nipype.utils.filemanip import split_filename
    from nipype.algorithms.misc import Gunzip
    from traits.trait_base import _Undefined

    if (in_file is None) or isinstance(in_file, _Undefined):
        return None

    if not isinstance(in_file, str):  # type(in_file) is list:
        return [unzip_nii(f) for f in in_file]

    _, base, ext = split_filename(in_file)

    # Not compressed
    if ext[-3:].lower() != ".gz":
        return in_file
    # Compressed
    gunzip = Gunzip(in_file=in_file)
    gunzip.run()
    return gunzip.aggregate_outputs().out_file


def fix_join(path, *paths):
    # This workaround is used in pipelines like DWIPreprocessingUsingT1
    # In the workflow.connect part, you can use some function that are used as string, causing an import error
    import os
    return os.path.join(path, *paths)


def check_bids_folder(bids_directory):
    import os
    from colorama import Fore
    from clinica.utils.exceptions import ClinicaBIDSError

    if not os.path.isdir(bids_directory):
        raise ClinicaBIDSError(
            "\n%s[Error] The BIDS directory you gave is not a folder.%s\n"
            "\n%sError explanations:%s\n"
            " - Clinica expected the following path to be a folder: %s%s%s\n"
            " - If you gave relative path, did you run Clinica on the good folder?" %
            (Fore.RED, Fore.RESET,
             Fore.YELLOW, Fore.RESET,
             Fore.BLUE, bids_directory, Fore.RESET)
        )


def caps_type_to_path(caps_type, caps_directory, participant_id, session_id):
    pipeline_to_path = {
        "t1-freesurfer": "%s/subjects/%s/%s/t1/freesurfer_cross_sectional/%s_%s" %
                         (caps_directory, participant_id, session_id, participant_id, session_id),
        "dwi-preprocessing": "%s/subjects/%s/%s/dwi/preprocessing" %
                             (caps_directory, participant_id, session_id),
    }
    dict_caps_type_to_path = {
        # t1-freesurfer
        "T1_FS_WM": "%s/mri/wm.seg.mgz" %
                    pipeline_to_path["t1-freesurfer"],
        "T1_FS_DESIKAN": "%s/mri/aparc+aseg.mgz" %
                         pipeline_to_path["t1-freesurfer"],
        "T1_FS_DESTRIEUX": "%s/mri/aparc.a2009s+aseg.mgz" %
                           pipeline_to_path["t1-freesurfer"],
        "T1_FS_BM": "%s/mri/brainmask.mgz" %
                    pipeline_to_path["t1-freesurfer"],
        "T1_FS_BE": "%s/mri/brain.mgz" %
                    pipeline_to_path["t1-freesurfer"],
        # dwi-preprocessing
        "DWI_PREPROC_BVAL": "%s/%s_%s*_dwi_space-(b0|T1w)_preproc.bval" %
                            (pipeline_to_path["dwi-preprocessing"], participant_id, session_id),
        "DWI_PREPROC_BVEC": "%s/%s_%s*_dwi_space-(b0|T1w)_preproc.bvec" %
                            (pipeline_to_path["dwi-preprocessing"], participant_id, session_id),
        "DWI_PREPROC_NII": "%s/%s_%s*_dwi_space-(b0|T1w)_preproc.(nii|nii.gz)" %
                           (pipeline_to_path["dwi-preprocessing"], participant_id, session_id),
        "DWI_PREPROC_BM": "%s/%s_%s*_dwi_space-(b0|T1w)_brainmask.(nii|nii.gz)" %
                          (pipeline_to_path["dwi-preprocessing"], participant_id, session_id),
    }
    return dict_caps_type_to_path[caps_type]


def caps_type_to_description(caps_type):
    dict_caps_type_to_description = {
        # t1-freesurfer
        "T1_FS_WM": "white matter segmentation",
        "T1_FS_DESIKAN": "Desikan parcellation",
        "T1_FS_DESTRIEUX": "Destrieux parcellation",
        "T1_FS_BM": "T1w brainmask",
        "T1_FS_BE": "brain extracted T1w",
        # dwi-preprocessing
        "DWI_PREPROC_BVAL": "preprocessed bval",
        "DWI_PREPROC_BVEC": "preprocessed bvec",
        "DWI_PREPROC_NII": "preprocessed DWI",
        "DWI_PREPROC_BM": "b0 brainmask",
    }
    return dict_caps_type_to_description[caps_type]


def check_input_caps_files(list_caps_files, caps_type, pipeline_name, caps_directory, participant_ids, session_ids):
    import re
    import collections
    from colorama import Fore

    subject_ids = [participant_ids[i] + '_' + session_ids[i] for i in range(len(participant_ids))]

    id_caps_files = [re.search(r'(sub-[a-zA-Z0-9]+)_(ses-[a-zA-Z0-9]+)', caps_file).group()
                     for caps_file in list_caps_files]

    id_duplicated_files = [item for item, count in collections.Counter(id_caps_files).items() if count > 1]

    id_missing_files = list(set(subject_ids) - set(id_caps_files))

    error_message = ""
    if id_missing_files:
        error_message += (
            "\n%s[Error] Clinica could not find %s file in CAPS directory for the following subject(s):%s\n" %
            (Fore.RED, caps_type_to_description(caps_type), Fore.RESET)
        )
        for id_missing_file in id_missing_files:
            participant_id = id_missing_file.split('_')[0]
            session_id = id_missing_file.split('_')[1]
            caps_path = caps_type_to_path(caps_type, caps_directory, participant_id, session_id)
            error_message += (
                "* %s%s%s. Expected path: %s\n" %
                (Fore.BLUE, id_missing_file.replace('_', '|'), Fore.RESET, caps_path)
            )
        error_message += (
            "%s* Error explanations:\n"
            " - Did you run the %s pipeline on the subject(s)?%s\n" %
            (Fore.YELLOW, pipeline_name, Fore.RESET)
        )

    if id_duplicated_files:
        error_message += (
            "\n%s[Error] Clinica found multiple %s files in BIDS directory for the following subject(s):%s" %
            (Fore.RED, caps_type_to_description(caps_type), Fore.RESET)
        )
        for id_duplicated_file in id_duplicated_files:
            error_message += (
                "\n* %s%s%s. Found files:\n%s" % (
                    Fore.BLUE, id_duplicated_file.replace('_', '|'), Fore.RESET,
                    '\n'.join('- ' + s for s in list_caps_files if id_duplicated_file in s)
                    # (Fore.BLUE + '|' + Fore.RESET).join(s for s in list_caps_files if id_duplicated_file in s)
                    ))
        error_message += (
            "\n%s* Error explanations:\n"
            " - Clinica expected to find a single file in CAPS directory for each subject.\n"
            " - Did you duplicate files in CAPS directory?%s\n" %
            (Fore.YELLOW, Fore.RESET)
        )
    return error_message
