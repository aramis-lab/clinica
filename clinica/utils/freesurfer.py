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


def get_secondary_stats(stats_filename, info_type):
    """Read the 'secondary' statistical info from .stats file

    Extract the information from .stats file that is commented out
    (lines starting with '# Measure' prefix) and does not appear in the
    table at the end of the document.

    Args:
        stats_filename (string): path to the .stats file
        info_type (string): 'volume', 'thickness' or 'area'

    Returns:
        secondary_stats_dict (dictionary of float): keys are regions of
            the brain, associated values are the corresponding
            volume/thickness/area depending on the input info type
    """
    # initialise structure containing the secondary statistical info
    secondary_stats_dict = dict()

    # currently no additional information is provided by .stats file for
    # the mean curvature
    if info_type != 'meancurv':
        # define how lines are supposed to end in the stats file, depending
        # on the type of information that is searched for
        endline_dict = dict()
        endline_dict['volume'] = 'mm^3'
        endline_dict['thickness'] = 'mm'
        endline_dict['area'] = 'mm^2'

        # define keywords that are supposed to appear in commented lines
        # containing statistical information
        info_keyword_dict = dict()
        info_keyword_dict['volume'] = ['volume', 'Volume']
        info_keyword_dict['area'] = ['area', 'Area']
        info_keyword_dict['thickness'] = ['thickness', 'Thickness']

        # read stats file line by line and only keep relevant lines
        with open(stats_filename, 'r') as stats_file:
            stats = stats_file.read()
        stats_line_list = stats.splitlines()
        for stats_line in stats_line_list:
            startswith_condition = stats_line.startswith('# Measure')
            endswith_condition = stats_line.endswith(endline_dict[info_type])
            if (startswith_condition) and (endswith_condition):
                stats_line_word_list = stats_line.replace(',', '').split()
                # sanity check: make sure any sensible variation of
                # 'volume', 'thickness' or 'area' appears inside the line
                if any(x in stats_line_word_list for x in info_keyword_dict[info_type]):
                    # add info
                    info_region = stats_line_word_list[2]
                    info_value = stats_line_word_list[-2]
                    secondary_stats_dict[info_region] = info_value

    return secondary_stats_dict


def generate_regional_measures(
        segmentation_path, subject_id, output_dir=None, longitudinal=False):
    """
    Read stats files located in
    <segmentation_path>/<subject_id>/stats/*.stats
    and generate TSV files in <segmentation_path>/regional_measures
    folder.

    Note: the .stats files contain both 1) a table with statistical
    information (e.g., structure volume) and 2) 'secondary' statistical
    information with all lines starting with the sentence '# Measure'.
    The .tsv files return the relevant statistical information from both
    sources.

    Args:
        segmentation_path (string): Path to the FreeSurfer segmentation.
        subject_id (string): Subject ID in the form sub-CLNC01_ses-M00
        output_dir (string): folder where the .tsv stats files will be
            stored. Will be [pass_segmentation]/regional_measures if no
            dir is provided by the user
        longitudinal (boolean): must be set to True if generating
            regional measures for a longitudinal processing (by default:
            set to False for a cross-sectional processing)
    """
    import os
    import errno
    import pandas
    from clinica.utils.stream import cprint
    from clinica.utils.freesurfer import write_tsv_file

    participant_id = subject_id.split('_')[0]
    session_id = subject_id.split('_')[1]

    # get location for recon-all stats files
    if not longitudinal:
        # cross-sectional processing
        stats_folder = os.path.join(segmentation_path, subject_id, 'stats')
    else:
        # longitudinal processing
        subses_long_id = '{0}.long.{1}'.format(subject_id, participant_id)
        stats_folder = os.path.join(segmentation_path, subses_long_id, 'stats')
    if not os.path.isdir(stats_folder):
        ioerror_msg = "Folder {0}/{1} does not contain FreeSurfer segmentation".format(
            segmentation_path, subject_id)
        raise IOError(ioerror_msg)

    # create output folder for the .tsv files
    # (define if not provided by user)
    if output_dir is None:
        output_dir = os.path.join(segmentation_path, 'regional_measures')
    try:
        os.makedirs(output_dir)
    except OSError as exception:
        # if dest_dir exists, go on, if its other error, raise
        if exception.errno != errno.EEXIST:
            raise

    # Generate TSV files for parcellation files
    #
    # Columns in ?h.BA.stats, ?h.aparc.stats or ?h.aparc.a2009s.stats file
    columns_parcellation = [
        'StructName', 'NumVert', 'SurfArea', 'GrayVol', 'ThickAvg', 'ThickStd',
        'MeanCurv', 'GausCurv', 'FoldInd', 'CurvInd'
    ]
    hemi_dict = {
        'left': 'lh',
        'right': 'rh'
    }
    atlas_dict = {
        'desikan': 'aparc',
        'destrieux': 'aparc.a2009s',
        'ba': 'BA_exvivo'
    }
    info_dict = {
        'volume': 'GrayVol',
        'thickness': 'ThickAvg',
        'area': 'SurfArea',
        'meancurv': 'MeanCurv'
    }
    for atlas in ('desikan', 'destrieux', 'ba'):
        stats_filename_dict = dict()
        df_dict = dict()
        # read both left and right .stats files
        for hemi in ('left', 'right'):
            stats_filename_dict[hemi] = os.path.join(
                stats_folder,
                '{0}.{1}.stats'.format(hemi_dict[hemi], atlas_dict[atlas]))
            df_dict[hemi] = pandas.read_csv(
                stats_filename_dict[hemi],
                names=columns_parcellation,
                comment='#', header=None, delimiter='\s+', dtype=str)
        # generate .tsv from 1) the table in .stats file and 2) the
        # secondary (commented out) information common to both 'left'
        # and 'right' .stats file
        for info in ('volume', 'thickness', 'area', 'meancurv'):
            # secondary information (common to 'left' and 'right')
            secondary_stats_dict = get_secondary_stats(
                stats_filename_dict['left'], info)
            # join primary and secondary information
            key_list = (
                list('lh_'+df_dict['left']['StructName'])
                + list('rh_'+df_dict['right']['StructName'])
                + list(secondary_stats_dict.keys()))
            col_name = info_dict[info]
            value_list = (
                list(df_dict['left'][col_name])
                + list(df_dict['right'][col_name])
                + list(secondary_stats_dict.values()))
            # Write .tsv
            write_tsv_file(
                os.path.join(
                    output_dir,
                    '{0}_parcellation-{1}_{2}.tsv'.format(
                        subject_id, atlas, info)),
                key_list,
                info,
                value_list)

    # Generate TSV files for segmentation files
    #
    # Columns in aseg.stats or wmparc.stats file
    columns_segmentation = [
        'Index', 'SegId', 'NVoxels', 'Volume_mm3', 'StructName', 'normMean',
        'normStdDev', 'normMin', 'normMax', 'normRange']

    # Parsing aseg.stats
    stats_filename = os.path.join(stats_folder, 'aseg.stats')
    df = pandas.read_csv(
        stats_filename,
        comment='#', header=None, delimiter='\s+', dtype=str,
        names=columns_segmentation)
    secondary_stats_dict = get_secondary_stats(stats_filename, 'volume')
    key_list = list(df['StructName'])+list(secondary_stats_dict.keys())
    value_list = list(df['Volume_mm3'])+list(secondary_stats_dict.values())
    write_tsv_file(
        os.path.join(output_dir, subject_id + '_segmentationVolumes.tsv'),
        key_list, 'volume', value_list)

    #  Parsing  wmparc.stats
    stats_filename = os.path.join(stats_folder, 'wmparc.stats')
    df = pandas.read_csv(
        stats_filename,
        comment='#', header=None, delimiter='\s+', dtype=str,
        names=columns_segmentation)
    secondary_stats_dict = get_secondary_stats(stats_filename, 'volume')
    key_list = list(df['StructName'])+list(secondary_stats_dict.keys())
    value_list = list(df['Volume_mm3'])+list(secondary_stats_dict.values())
    write_tsv_file(
        os.path.join(output_dir, subject_id + '_parcellation-wm_volume.tsv'),
        key_list, 'volume', value_list)


def write_tsv_file(out_filename, name_list, scalar_name, scalar_list):
    """Write a .tsv file with list of keys and values

    Args:
        out_filename (string): name of the .tsv file
        name_list (list of string): list of keys
        scalar_name (string): 'volume', 'thickness', 'area' or
            'meanCurv'. Not used for now. Might be used as part of a
            pandas data frame
        scalar_list (list of float): list of values corresponding to the keys
    """
    import pandas
    import warnings

    try:
        data = pandas.DataFrame({
            'label_name': name_list,
            'label_value': scalar_list})
        data.to_csv(out_filename, sep='\t', index=False, encoding='utf-8')
    except Exception as exception:
        warnings.warn("Impossible to save {0} file".format(out_filename))
        raise exception

    return out_filename
