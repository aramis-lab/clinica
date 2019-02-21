# coding: utf8


"""This module contains FreeSurfer utilities."""


def freesurfer_volume_to_native_volume(
        freesurfer_volume,
        native_volume,
        name_output_volume=None):
    """
    Convert FreeSurfer volume in native space.

    This function converts any volume in FreeSurfer's conformed space
    (1x1x1mm voxel size, 256x256x256 dimension) into a volume in native space.

    For further details:
    https://surfer.nmr.mgh.harvard.edu/fswiki/FsAnat-to-NativeAnat

    Args:
        freesurfer_volume (str): Volume in FreeSurfer's conformed space
            (e.g. aparc+aseg.mgz containing the Desikan parcellation)
        native_volume (str): Volume in native space (You should choose
            ${SUBJECTS_DIR}/subject_id/mri/rawavg.mgz).
        name_output_volume (Optional[str]): Name of the output matrix
            (default=volume_in_native_space.nii.gz).

    Returns:
        out_volume (str): volume in native space (the file is saved here:
            ${SUBJECTS_DIR}/subject_id/native_space/label.nii)

    Example:
        >>> from clinica.utils.freesurfer import freesurfer_volume_to_native_volume
        >>> freesurfer_volume_to_native_volume(bert/mri/rawavg.mgz, bert/mri/aparc+aseg.mgz, 'aparc-in-native-space.nii')
    """
    import os
    import os.path as op

    assert(op.isfile(freesurfer_volume))
    assert(op.isfile(native_volume))

    if name_output_volume is None:
        out_volume = op.abspath('volume_in_native_space.nii.gz')
    else:
        out_volume = op.abspath(name_output_volume)

    cmd = 'mri_vol2vol --regheader --no-save-reg --mov %s --targ %s --o %s' \
          % (freesurfer_volume, native_volume, out_volume)
    os.system(cmd)

    return out_volume


def fs_caps2reconall(caps_dir, dest_dir, subjects_visits_tsv):
    """
    This function transfers CAPS recon-all output structure to
    standard FreeSurfer recon-all output structure.

    Args:
        caps_dir: CAPS directory containing the outputs in CAPS hierarchy
        dest_dir: the destination folder containing the FreeSurfer output structure
        subjects_visits_tsv: tsv files containing the subjects that you want
            to convert
    """
    import os
    import csv
    from shutil import copytree
    from clinica.utils.stream import cprint

    subject_list = []
    session_list = []
    with open(subjects_visits_tsv, 'rb') as tsvin:
        tsv_reader = csv.reader(tsvin, delimiter='\t')

        for row in tsv_reader:
            if row[0] == 'participant_id':
                continue
            else:
                subject_list.append(row[0])
                session_list.append(row[1])

    output_path = os.path.expanduser(caps_dir)  # change the relative path to be absolute path
    caps_dir = os.path.join(output_path, 'subjects')

    for i in range(len(subject_list)):
        if os.path.isdir(os.path.join(dest_dir, subject_list[i] + '_' + session_list[i])):
            cprint("This subject: %s for FreeSurfer exits already!" % subject_list[i])
        else:
            cprint("Convert subject: %s from CAPS to FreeSurfer output structure" % subject_list[i])
            copytree(os.path.join(caps_dir, subject_list[i], session_list[i], 't1/freesurfer_cross_sectional', subject_list[i] + '_' + session_list[i]), os.path.join(dest_dir, subject_list[i] + '_' + session_list[i]))
            cprint("--------------Finish this subject!-----------------------")


def write_volumetric_per_subject(caps_dir, subjects_visits_tsv):
    """
        This func is to write the volumetric measurement after recon-all
        pipelines for each subjects in the subjects_visits_tsv

    Args: caps_dir: CAPS directory subjects_visits_tsv: tsv contains all the
    particiapnt_id and session_id

    Returns:

    """
    import nipype.pipeline.engine as pe
    from nipype.interfaces.utility import Function
    import pandas as pd
    from clinica.pipelines.t1_freesurfer_cross_sectional.t1_freesurfer_cross_sectional_utils import write_statistics_per_subject

    # get the list for subject_ids
    subjects_visits = pd.io.parsers.read_csv(subjects_visits_tsv, sep='\t')
    if (list(subjects_visits.columns.values)[0] != 'participant_id') and (list(subjects_visits.columns.values)[1] != 'session_id'):
        raise Exception('Subjects and visits file is not in the correct format.')
    subject_list = list(subjects_visits.participant_id)
    session_list = list(subjects_visits.session_id)
    subject_id = list(subject_list[i] + '_' + session_list[i] for i in range(len(subject_list)))

    fs_tsv_subject = pe.MapNode(name='volumetric_summary_node',
                                iterfield=['subject_id'],
                                interface=Function(
                                    input_names=['subject_id', 'output_dir'],
                                    output_names=[],
                                    function=write_statistics_per_subject,
                                    imports=['import os', 'import errno']))
    fs_tsv_subject.inputs.subject_id = subject_id
    fs_tsv_subject.inputs.output_dir = caps_dir

    return fs_tsv_subject


def write_reconall_log_summary(caps_dir, subjects_visits_tsv):
    """
        This func is to write the recon_all.log summary for all the subjects,
        the first step quality check

    Args: caps_dir: CAPS directory subjects_visits_tsv: tsv contains all the
    particiapnt_id and session_id

    Returns:

    """
    import nipype.pipeline.engine as pe
    from nipype.interfaces.utility import Function
    import pandas as pd
    from clinica.pipelines.t1_freesurfer_cross_sectional.t1_freesurfer_cross_sectional_utils import log_summary

    # get the list for subject_ids
    subjects_visits = pd.io.parsers.read_csv(subjects_visits_tsv, sep='\t')
    if ((list(subjects_visits.columns.values)[0] != 'participant_id')
            and (list(subjects_visits.columns.values)[1] != 'session_id')):
        raise Exception('Subjects and visits file is not in the correct format.')
    subject_list = list(subjects_visits.participant_id)
    session_list = list(subjects_visits.session_id)
    subject_id = list(subject_list[i] + '_' + session_list[i] for i in range(len(subject_list)))

    lognode = pe.Node(name='lognode',
                      interface=Function(
                          input_names=['subject_list', 'session_list', 'subject_id', 'output_dir'],
                          output_names=[],
                          function=log_summary))
    lognode.inputs.subject_list = subject_list
    lognode.inputs.session_list = session_list
    lognode.inputs.subject_id = subject_id
    lognode.inputs.output_dir = caps_dir

    return lognode


def generate_regional_measures(path_segmentation, subject_id, output_dir=None):
    """
    Read stats files located in <path_segmentation>/<subject_id>/stats/*.stats
    and generate TSV files in <path_segmentation>/regional_measures folder.

    Args:
        path_segmentation: Path to the FreeSurfer segmentation.
        subject_id: Subject ID in the form sub-CLNC01_ses-M00
        output_dir: CAPS directory
    """
    import os
    import errno
    import pandas
    from clinica.utils.stream import cprint
    from clinica.utils.freesurfer import write_tsv_file

    participant_id = subject_id.split('_')[0]
    session_id = subject_id.split('_')[1]

    stats_folder = os.path.join(path_segmentation, subject_id, 'stats')
    if not os.path.isdir(stats_folder):
        raise IOError("Folder %s/%s does not contain FreeSurfer segmentation" % (path_segmentation, subject_id))

    if output_dir is None:
        output_dir = os.path.join(path_segmentation, 'regional_measures')

    try:
        os.makedirs(output_dir)
    except OSError as exception:
        if exception.errno != errno.EEXIST:  # if dest_dir exists, go on, if its other error, raise
            raise

    #
    # Generate TSV files for parcellation files
    #
    # Columns in ?h.BA.stats, ?h.aparc.stats or ?h.aparc.a2009s.stats file
    columns_parcellation = [
        'StructName', 'NumVert', 'SurfArea', 'GrayVol', 'ThickAvg', 'ThickStd',
        'MeanCurv', 'GausCurv', 'FoldInd', 'CurvInd'
    ]
    dict_hemi = {
        "left": "lh",
        "right": "rh"
    }
    dict_atlas = {
        "desikan": "aparc",
        "destrieux": "aparc.a2009s",
        "ba": "BA_exvivo"
    }
    for atlas in ('desikan', 'destrieux', 'ba'):
        for hemi in ('left', 'right'):
            df = pandas.read_csv(
                os.path.join(stats_folder, dict_hemi[hemi]+'.' + dict_atlas[atlas]+'.stats'),
                names=columns_parcellation,
                comment='#', header=None, delimiter='\s+', dtype=str)
            write_tsv_file(os.path.join(output_dir, subject_id+'_hemi-'+hemi+'_parcellation-'+atlas+'_volume.tsv'),
                           df['StructName'], 'volume', df['GrayVol'])
            write_tsv_file(os.path.join(output_dir, subject_id+'_hemi-'+hemi+'_parcellation-'+atlas+'_thickness.tsv'),
                           df['StructName'], 'thickness', df['ThickAvg'])
            write_tsv_file(os.path.join(output_dir, subject_id+'_hemi-'+hemi+'_parcellation-'+atlas+'_area.tsv'),
                           df['StructName'], 'area', df['SurfArea'])
            write_tsv_file(os.path.join(output_dir, subject_id+'_hemi-'+hemi+'_parcellation-'+atlas+'_meancurv.tsv'),
                           df['StructName'], 'meancurv', df['MeanCurv'])

    #
    # Generate TSV files for segmentation files
    #
    # Columns in aseg.stats or wmparc.stats file
    columns_segmentation = [
        'Index', 'SegId', 'NVoxels', 'Volume_mm3', 'StructName', 'normMean',
        'normStdDev', 'normMin', 'normMax', 'normRange']

    # Parsing aseg.stats
    df = pandas.read_csv(os.path.join(stats_folder, 'aseg.stats'),
                         comment='#', header=None, delimiter='\s+', dtype=str,
                         names=columns_segmentation)
    write_tsv_file(os.path.join(output_dir, subject_id + '_segmentationVolumes.tsv'),
                   df['StructName'], 'volume', df['Volume_mm3'])

    #  Parsing  wmparc.stats
    df = pandas.read_csv(os.path.join(stats_folder, 'wmparc.stats'),
                         comment='#', header=None, delimiter='\s+', dtype=str,
                         names=columns_segmentation)
    write_tsv_file(os.path.join(output_dir, subject_id + '_parcellation-wm_volume.tsv'),
                   df['StructName'], 'volume', df['Volume_mm3'])


def write_tsv_file(out_file, list_names, name_scalar, list_scalars):
    # import pandas
    import csv
    from clinica.utils.stream import cprint

    try:
        with open(out_file, 'w', encoding='utf8', newline='') as tsv_file:
            tsv_writer = csv.writer(tsv_file, delimiter='\t',
                                    lineterminator='\n')
            tsv_writer.writerow(list_names)
            tsv_writer.writerow(list_scalars)
        # TODO: Save data in lines instead of columns for Clinica 0.3
        # data = pandas.DataFrame({
        #     'label_name': list_names,
        #     name_scalar: list_scalars
        # })
        # data.to_csv(out_file, sep='\t', index=True, encoding='utf-8')
    except Exception as e:
        cprint("Impossible to save %s file" % out_file)
        raise e

    return out_file
