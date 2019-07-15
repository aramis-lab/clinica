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
    # This workaround is used in pipelines like DWIPreporcessingUsingT1
    # In the workflow.connect part, you can use some function that are used as string, causing an import error
    import os
    return os.path.join(path, *paths)


def check_input_bids_file(list_bids_files, bids_type, bids_directory, participant_id, session_id):
    from clinica.utils.exceptions import ClinicaBIDSError
    from colorama import Fore

    bids_type_to_description = {
        # anat/
        "T1W_NII": "T1 weighted image",
        # dwi/
        "DWI_BVAL": "bval file",
        "DWI_BVEC": "bvec file",
        "DWI_NII": "diffusion weighted image",
        "DWI_JSON": "DWI JSON file",
        # fmap/
        "FMAP_MAGNITUDE1_NII": "magnitude (1st) image",
        "FMAP_MAGNITUDE1_JSON": "magnitude (1st) JSON file",
        "FMAP_MAGNITUDE2_NII": "magnitude (2st) image",
        "FMAP_MAGNITUDE2_JSON": "magnitude (2st) JSON file",
        "FMAP_PHASEDIFF_NII": "phase difference image",
        "FMAP_PHASEDIFF_JSON": "phase difference JSON file",
        "FMAP_PHASE1_NII": "phase (1st) image",
        "FMAP_PHASE1_JSON": "phase (1st) JSON file",
        "FMAP_PHASE2_NII": "phase (2nd) image",
        "FMAP_PHASE2_JSON": "phase (2nd) JSON file",
        "FMAP_MAGNITUDE_NII": "magnitude image",
        "FMAP_MAGNITUDE_JSON": "magnitude JSON file",
        "FMAP_NII": "fieldmap image",
        "FMAP_JSON": "fieldmap JSON file",
    }

    bids_type_to_path = {
        # anat/
        "T1W_NII": "%s/%s/anat/%s_%s*_T1w.(nii|nii.gz)" % (participant_id, session_id, participant_id, session_id),
        # dwi/
        "DWI_BVAL": "%s/%s/dwi/%s_%s*_dwi.bval" % (participant_id, session_id, participant_id, session_id),
        "DWI_BVEC": "%s/%s/dwi/%s_%s*_dwi.bvec" % (participant_id, session_id, participant_id, session_id),
        "DWI_NII": "%s/%s/dwi/%s_%s*_dwi.(nii|nii.gz)" % (participant_id, session_id, participant_id, session_id),
        "DWI_JSON": "%s/%s/dwi/%s_%s*_dwi.json" % (participant_id, session_id, participant_id, session_id),
        # fmap/
        "FMAP_MAGNITUDE1_NII": "%s/%s/dwi/%s_%s*_magnitude1.(nii|nii.gz)" % (participant_id, session_id, participant_id, session_id),
        "FMAP_MAGNITUDE1_JSON": "%s/%s/dwi/%s_%s*_magnitude1.json" % (participant_id, session_id, participant_id, session_id),
        "FMAP_MAGNITUDE2_NII": "%s/%s/dwi/%s_%s*_magnitude2.(nii|nii.gz)" % (participant_id, session_id, participant_id, session_id),
        "FMAP_MAGNITUDE2_JSON": "%s/%s/dwi/%s_%s*_magnitude2.json" % (participant_id, session_id, participant_id, session_id),
        "FMAP_PHASEDIFF_NII": "%s/%s/dwi/%s_%s*_phasediff.(nii|nii.gz)" % (participant_id, session_id, participant_id, session_id),
        "FMAP_PHASEDIFF_JSON": "%s/%s/dwi/%s_%s*_phasediff.json" % (participant_id, session_id, participant_id, session_id),
        "FMAP_PHASE1_NII": "%s/%s/dwi/%s_%s*_phase1.(nii|nii.gz)" % (participant_id, session_id, participant_id, session_id),
        "FMAP_PHASE1_JSON": "%s/%s/dwi/%s_%s*_phase1.json" % (participant_id, session_id, participant_id, session_id),
        "FMAP_PHASE2_NII": "%s/%s/dwi/%s_%s*_phase2.(nii|nii.gz)" % (participant_id, session_id, participant_id, session_id),
        "FMAP_PHASE2_JSON": "%s/%s/dwi/%s_%s*_phase2.json" % (participant_id, session_id, participant_id, session_id),
        "FMAP_MAGNITUDE_NII": "%s/%s/dwi/%s_%s*_magnitude.(nii|nii.gz)" % (participant_id, session_id, participant_id, session_id),
        "FMAP_MAGNITUDE_JSON": "%s/%s/dwi/%s_%s*_magnitude.json" % (participant_id, session_id, participant_id, session_id),
        "FMAP_NII": "%s/%s/dwi/%s_%s*_fieldmap.(nii|nii.gz)" % (participant_id, session_id, participant_id, session_id),
        "FMAP_JSON": "%s/%s/dwi/%s_%s*_fieldmap.json" % (participant_id, session_id, participant_id, session_id),
    }

    if len(list_bids_files) == 0:
        raise ClinicaBIDSError(
            "\n%s[Error] Clinica could not find %s file in BIDS directory for participant %s at session %s%s.\n"
            "\n%sError explanations:%s\n"
            " - Clinica expected to find the file at the following path: %s%s/%s%s\n"
            " - Did the subject have the acquisition?" %
            (Fore.RED, bids_type_to_description[bids_type], participant_id[4:], session_id[4:], Fore.RESET,
             Fore.YELLOW, Fore.RESET,
             Fore.BLUE, bids_directory, bids_type_to_path[bids_type], Fore.RESET)
        )
    elif len(list_bids_files) > 1:
        raise ClinicaBIDSError(
            "\n%s[Error] Clinica found %s %s files in CAPS directory for participant %s at session %s.%s\n"
            "\n%sError explanations:%s\n"
            " - Clinica expected to find a single file in BIDS directory. Found files:%s%s%s\n"
            " - If you have different runs and/or acquisitions on the same folder, Clinica can not handle this situation." %
            (Fore.RED, len(list_bids_files), bids_type_to_description[bids_type], participant_id[4:], session_id[4:], Fore.RESET,
             Fore.YELLOW, Fore.RESET,
             Fore.BLUE, list_bids_files, Fore.RESET)
        )
    if len(list_bids_files) == 0:
        raise ClinicaBIDSError(
            "Missing %s in BIDS dataset for participant %s%s%s at session %s%s%s.\n%sExpected path:\n%s%s/%s" %
            (bids_type_to_description[bids_type],
             Fore.BLUE, participant_id[4:], Fore.RESET,
             Fore.BLUE, session_id[4:], Fore.RESET,
             Fore.YELLOW, Fore.RESET, bids_directory, bids_type_to_path[bids_type])
        )
    elif len(list_bids_files) > 1:
        raise ClinicaBIDSError(
            "Found %s %ss in BIDS dataset for %s at session %s: you should only have one file." %
            (len(list_bids_files), bids_type_to_description[bids_type], participant_id, session_id)
        )


def check_input_caps_file(list_caps_files, caps_type, pipeline_name, caps_directory, participant_id, session_id):
    from clinica.utils.exceptions import ClinicaCAPSError
    from colorama import Fore
    caps_type_to_description = {
        # t1-freesurfer
        "T1_FS_WM": "white matter segmentation",
        "T1_FS_DESIKAN": "Desikan parcellation",
        "T1_FS_DESTRIEUX": "Destrieux parcellation",
        "T1_FS_BM": "T1w brainmask",
        # dwi-preprocessing
        "DWI_PREPROC_BVAL": "preprocessed bval",
        "DWI_PREPROC_BVEC": "preprocessed bvec",
        "DWI_PREPROC_NII": "preprocessed DWI",
        "DWI_PREPROC_BM": "b0 brainmask",
    }
    pipeline_to_path = {
        "t1-freesurfer": "subjects/%s/%s/t1/freesurfer_cross_sectional/%s_%s" % (participant_id, session_id, participant_id, session_id),
        "dwi-preprocessing": "subjects/%s/%s/dwi/preprocessing" % (participant_id, session_id),
    }
    caps_type_to_path = {
        # t1-freesurfer
        "T1_FS_WM": "%s/mri/wm.seg.mgz" % pipeline_to_path["t1-freesurfer"],
        "T1_FS_DESIKAN": "%s/mri/aparc+aseg.mgz" % pipeline_to_path["t1-freesurfer"],
        "T1_FS_DESTRIEUX": "%s/mri/aparc.a2009s+aseg.mgz" % pipeline_to_path["t1-freesurfer"],
        "T1_FS_BM": "%s/mri/brainmask.mgz" % pipeline_to_path["t1-freesurfer"],
        # dwi-preprocessing
        "DWI_PREPROC_BVAL": "%s/%s_%s*_preproc.bval" % (pipeline_to_path["dwi-preprocessing"], participant_id, session_id),
        "DWI_PREPROC_BVEC": "%s/%s_%s*_preproc.bvec" % (pipeline_to_path["dwi-preprocessing"], participant_id, session_id),
        "DWI_PREPROC_NII": "%s/%s_%s*_preproc.(nii|nii.gz)" % (pipeline_to_path["dwi-preprocessing"], participant_id, session_id),
        "DWI_PREPROC_BM": "%s/%s_%s*_brainmask.(nii|nii.gz)" % (pipeline_to_path["dwi-preprocessing"], participant_id, session_id),
    }

    if len(list_caps_files) == 0:
        raise ClinicaCAPSError(
            "\n%s[Error] Clinica could not find %s file in CAPS directory for participant %s at session %s.%s\n"
            "\n%sError explanations:%s\n"
            " - Clinica expected to find the file at the following path: %s%s/%s%s\n"
            " - Did you run the %s%s%s pipeline on this image?" %
            (Fore.RED, caps_type_to_description[caps_type], participant_id[4:], session_id[4:], Fore.RESET,
             Fore.YELLOW, Fore.RESET,
             Fore.BLUE, caps_directory, caps_type_to_path[caps_type], Fore.RESET,
             Fore.BLUE, pipeline_name, Fore.RESET)
        )
    elif len(list_caps_files) > 1:
        raise ClinicaCAPSError(
            "\n%s[Error] Clinica found %s %s files in CAPS directory for participant %s at session %s%s.\n"
            "\n%sError explanations:%s\n"
            " - Clinica expected to find a single file in CAPS directory. Found files:%s%s%s\n"
            " - Did you duplicate files in CAPS directory: %s%s/%s%s?" %
            (Fore.RED, len(list_caps_files), caps_type_to_description[caps_type], participant_id[4:], session_id[4:], Fore.RESET,
             Fore.YELLOW, Fore.RESET,
             Fore.BLUE, list_caps_files, Fore.RESET,
             Fore.BLUE, caps_directory, pipeline_to_path[pipeline_name], Fore.RESET)
        )
