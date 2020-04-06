# coding: utf8

"""
Describe files to grab, to use with inputs.clinica_file_reader() and inputs.clinica_group_reader()

"""

""" T1w """

# BIDS

T1W_NII = {'pattern': 'sub-*_ses-*_t1w.nii*',
           'description': 'T1w MRI'}

# T1-FreeSurfer

T1_FS_WM = {'pattern': 't1/freesurfer_cross_sectional/sub-*_ses-*/mri/wm.seg.mgz',
            'description': 'segmentation of white matter (mri/wm.seg.mgz).',
            'needed_pipeline': 't1-freesurfer'}

T1_FS_BRAIN = {'pattern': 't1/freesurfer_cross_sectional/sub-*_ses-*/mri/brain.mgz',
               'description': ' extracted brain from T1w MRI (mri/brain.mgz).',
               'needed_pipeline': 't1-freesurfer'}

T1_FS_ORIG_NU = {'pattern': 't1/freesurfer_cross_sectional/sub-*_ses-*/mri/orig_nu.mgz',
                 'description': 'intensity normalized volume generated after correction for'
                                ' non-uniformity in FreeSurfer (mri/orig_nu.mgz).',
                 'needed_pipeline': 't1-freesurfer'}

T1_FS_LONG_ORIG_NU = {'pattern': 't1/long-*/freesurfer_longitudinal/sub-*_ses-*.long.sub-*_*/mri/orig_nu.mgz',
                      'description': 'intensity normalized volume generated after correction for non-uniformity in FreeSurfer (orig_nu.mgz) in longitudinal',
                      'needed_pipeline': 't1-freesurfer and t1-freesurfer longitudinal'}

T1_FS_WM_SURF_R = {'pattern': 't1/freesurfer_cross_sectional/sub-*_ses-*/surf/rh.white',
                   'description': 'right white matter/gray matter border surface (rh.white).',
                   'needed_pipeline': 't1-freesurfer'}

T1_FS_LONG_SURF_R = {'pattern': 't1/long-*/freesurfer_longitudinal/sub-*_ses-*.long.sub-*_*/surf/rh.white',
                     'description': 'right white matter/gray matter border surface (rh.white) generated with t1-freesurfer-longitudinal.',
                     'needed_pipeline': 't1-freesurfer and t1-freesurfer longitudinal'}

T1_FS_LONG_SURF_L = {'pattern': 't1/long-*/freesurfer_longitudinal/sub-*_ses-*.long.sub-*_*/surf/lh.white',
                     'description': 'left white matter/gray matter border surface (lh.white) generated with t1-freesurfer-longitudinal.',
                     'needed_pipeline': 't1-freesurfer and t1-freesurfer longitudinal'}

T1_FS_WM_SURF_L = {'pattern': 't1/freesurfer_cross_sectional/sub-*_ses-*/surf/lh.white',
                   'description': 'left white matter/gray matter border surface (lh.white).',
                   'needed_pipeline': 't1-freesurfer'}

T1_FS_DESTRIEUX = {'pattern': 't1/freesurfer_cross_sectional/sub-*_ses-*/mri/aparc.a2009s+aseg.mgz',
                   'description': 'Destrieux-based segmentation (mri/aparc.a2009s+aseg.mgz).',
                   'needed_pipeline': 't1-freesurfer'}

T1_FS_DESTRIEUX_PARC_L = {'pattern': 't1/freesurfer_cross_sectional/sub-*_ses-*/label/lh.aparc.a2009s.annot',
                          'description': 'left hemisphere surface-based Destrieux parcellation (label/lh.aparc.a2009s.annot).',
                          'needed_pipeline': 't1-freesurfer'}

T1_FS_LONG_DESTRIEUX_PARC_L = {'pattern': 't1/long-*/freesurfer_longitudinal/sub-*_ses-*.long.sub-*_*/label/lh.aparc.a2009s.annot',
                               'description': 'left hemisphere surface-based Destrieux parcellation (label/lh.aparc.a2009s.annot) generated with t1-freesurfer-longitudinal.',
                               'needed_pipeline': 't1-freesurfer and t1-freesurfer longitudinal'}

T1_FS_LONG_DESTRIEUX_PARC_R = {'pattern': 't1/long-*/freesurfer_longitudinal/sub-*_ses-*.long.sub-*_*/label/rh.aparc.a2009s.annot',
                               'description': 'right hemisphere surface-based Destrieux parcellation (label/rh.aparc.a2009s.annot) generated with t1-freesurfer-longitudinal.',
                               'needed_pipeline': 't1-freesurfer and t1-freesurfer longitudinal'}

T1_FS_DESTRIEUX_PARC_R = {'pattern': 't1/freesurfer_cross_sectional/sub-*_ses-*/label/rh.aparc.a2009s.annot',
                          'description': 'right hemisphere surface-based Destrieux parcellation (label/rh.aparc.a2009s.annot).',
                          'needed_pipeline': 't1-freesurfer'}

T1_FS_DESIKAN = {'pattern': 't1/freesurfer_cross_sectional/sub-*_ses-*/mri/aparc+aseg.mgz',
                 'description': 'Desikan-based segmentation (mri/aparc.a2009s+aseg.mgz).',
                 'needed_pipeline': 't1-freesurfer'}

T1_FS_DESIKAN_PARC_L = {'pattern': 't1/freesurfer_cross_sectional/sub-*_ses-*/label/lh.aparc.annot',
                        'description': 'left hemisphere surface-based Desikan parcellation (label/lh.aparc.annot).',
                        'needed_pipeline': 't1-freesurfer'}

T1_FS_DESIKAN_PARC_R = {'pattern': 't1/freesurfer_cross_sectional/sub-*_ses-*/label/rh.aparc.annot',
                        'description': 'right hemisphere surface-based Desikan parcellation (label/rh.aparc.annot).',
                        'needed_pipeline': 't1-freesurfer'}

# T1-FreeSurfer-Template
T1_FS_T_DESTRIEUX = {'pattern': 'freesurfer_unbiased_template/sub-*_long-*/mri/aparc.a2009s+aseg.mgz',
                     'description': 'Destrieux-based segmentation (mri/aparc.a2009s+aseg.mgz) from unbiased template.',
                     'needed_pipeline': 't1-freesurfer-longitudinal or t1-freesurfer-template'}

# T1-FreeSurfer-Longitudinal-Correction
T1_FS_LONG_DESIKAN_PARC_L = {'pattern': 't1/long-*/freesurfer_longitudinal/sub-*_ses-*.long.sub-*_*/label/lh.aparc.annot',
                             'description': 'left hemisphere surface-based Desikan parcellation (label/lh.aparc.annot) generated with t1-freesurfer-longitudinal.',
                             'needed_pipeline': 't1-freesurfer and t1-freesurfer-longitudinal'}

T1_FS_LONG_DESIKAN_PARC_R = {'pattern': 't1/long-*/freesurfer_longitudinal/sub-*_ses-*.long.sub-*_*/label/rh.aparc.annot',
                             'description': 'right hemisphere surface-based Desikan parcellation (label/rh.aparc.annot) generated with t1-freesurfer-longitudinal.',
                             'needed_pipeline': 't1-freesurfer and t1-freesurfer-longitudinal'}


# T1-Volume
def t1_volume_native_tpm(tissue_number):
    from .spm import INDEX_TISSUE_MAP
    import os
    information = {
        'pattern': os.path.join('t1', 'spm', 'segmentation', 'native_space',
                                '*_*_T1w_segm-' + INDEX_TISSUE_MAP[tissue_number] + '_probability.nii*'),
        'description': 'Tissue probability map ' + INDEX_TISSUE_MAP[tissue_number] + ' in native space',
        'needed_pipeline': 't1-volume-tissue-segmentation'
    }
    return information


def t1_volume_dartel_input_tissue(tissue_number):
    from .spm import INDEX_TISSUE_MAP
    import os
    information = {
        'pattern': os.path.join('t1', 'spm', 'segmentation', 'dartel_input',
                                '*_*_T1w_segm-' + INDEX_TISSUE_MAP[tissue_number] + '_dartelinput.nii*'),
        'description': 'Dartel input for tissue probability map ' + INDEX_TISSUE_MAP[tissue_number] + ' from T1w MRI',
        'needed_pipeline': 't1-volume-tissue-segmentation'
    }
    return information


def t1_volume_native_tpm_in_mni(tissue_number, modulation):
    from .spm import INDEX_TISSUE_MAP
    import os
    if modulation:
        pattern_modulation = 'on'
        description_modulation = "with"
    else:
        pattern_modulation = 'off'
        description_modulation = "without"
    information = {
        'pattern': os.path.join('t1', 'spm', 'segmentation', 'normalized_space',
                                '*_*_T1w_segm-' + INDEX_TISSUE_MAP[tissue_number] + '_space-Ixi549Space' +
                                '_modulated-' + pattern_modulation + '_probability.nii*'),
        'description': 'Tissue probability map ' + INDEX_TISSUE_MAP[tissue_number] + ' based on native MRI ' +
                       'in MNI space (Ixi549) ' + description_modulation + ' modulation.',
        'needed_pipeline': 't1-volume-tissue-segmentation'}
    return information


def t1_volume_template_tpm_in_mni(group_label, tissue_number, modulation):
    from .spm import INDEX_TISSUE_MAP
    import os
    if modulation:
        pattern_modulation = 'on'
        description_modulation = "with"
    else:
        pattern_modulation = 'off'
        description_modulation = "without"
    information = {
        'pattern': os.path.join('t1', 'spm', 'dartel', 'group-' + group_label,
                                '*_T1w_segm-' + INDEX_TISSUE_MAP[tissue_number] + '_space-Ixi549Space_modulated-' +
                                pattern_modulation + '_probability.nii*'),
        'description': 'Tissue probability map ' + INDEX_TISSUE_MAP[tissue_number] + ' based on ' + group_label +
                       ' template in MNI space (Ixi549) ' + description_modulation + ' modulation.',
        'needed_pipeline': 't1-volume'
    }
    return information


def t1_volume_deformation_to_template(group_label):
    import os
    information = {
        'pattern': os.path.join('t1', 'spm', 'dartel', 'group-' + group_label,
                                'sub-*_ses-*_T1w_target-' + group_label + '_transformation-forward_deformation.nii*'),
        'description': 'Deformation from native space to group template ' + group_label + ' space.',
        'needed_pipeline': 't1-volume-create-dartel'
    }
    return information


def t1_volume_i_th_iteration_group_template(group_label, i):
    import os
    information = {
        'pattern': os.path.join('group-' + group_label, 't1',
                                'group-' + group_label + '_iteration-' + str(i) + '_template.nii*'),
        'description': 'Iteration #' + str(i) + ' of Dartel template ' + group_label,
        'needed_pipeline': 't1-volume or t1-volume-create-dartel'
    }
    return information


def t1_volume_final_group_template(group_label):
    import os
    information = {
        'pattern': os.path.join('group-' + group_label, 't1',
                                'group-' + group_label + '_template.nii*'),
        'description': 'T1w template file of group ' + group_label,
        'needed_pipeline': 't1-volume or t1-volume-create-dartel'
    }
    return information


""" DWI """

# BIDS

DWI_NII = {'pattern': 'dwi/sub-*_ses-*_dwi.nii*',
           'description': 'DWI NIfTI'}

DWI_JSON = {'pattern': 'dwi/sub-*_ses-*_dwi.json',
            'description': 'DWI JSON file'}

DWI_BVAL = {'pattern': 'dwi/sub-*_ses-*_dwi.bval',
            'description': 'bval files'}

DWI_BVEC = {'pattern': 'dwi/*_dwi.bvec',
            'description': 'bvec files'}

FMAP_PHASEDIFF_JSON = {'pattern': 'fmap/sub-*_ses-*_phasediff.json',
                       'description': 'phasediff JSON file'}

FMAP_PHASEDIFF_NII = {'pattern': 'fmap/sub-*_ses-*_phasediff.nii*',
                      'description': 'phasediff NIfTI volume'}

FMAP_MAGNITUDE1_NII = {'pattern': 'fmap/sub-*_ses-*_magnitude1.nii*',
                       'description': 'magnitude1 file'}

# CAPS

DWI_PREPROC_NII = {'pattern': 'dwi/preprocessing/sub-*_ses-*_dwi_space-*_preproc.nii*',
                   'description': 'preprocessed DWI',
                   'needed_pipeline': 'dwi-preprocessing-using-t1 or dwi-preprocessing-using-fieldmap'}

DWI_PREPROC_BRAINMASK = {'pattern': 'dwi/preprocessing/sub-*_ses-*_dwi_space-*_brainmask.nii*',
                         'description': 'b0 brainmask',
                         'needed_pipeline': 'dwi-preprocessing-using-t1 or dwi-preprocessing-using-fieldmap'}

DWI_PREPROC_BVEC = {'pattern': 'dwi/preprocessing/sub-*_ses-*_dwi_space-*_preproc.bvec',
                    'description': 'preprocessed bvec',
                    'needed_pipeline': 'dwi-preprocessing-using-t1 or dwi-preprocessing-using-fieldmap'}

DWI_PREPROC_BVAL = {'pattern': 'dwi/preprocessing/*_dwi_space-*_preproc.bval',
                    'description': 'preprocessed bval',
                    'needed_pipeline': 'dwi-preprocessing-using-t1 or dwi-preprocessing-using-fieldmap'}

""" fMRI """

# BIDS
FMRI_BOLD_NII = {'pattern': 'func/sub-*_ses-*_bold.nii*',
                 'description': 'blood-oxygen-level-dependent (BOLD) image'}

FMRI_BOLD_JSON = {'pattern': 'func/sub-*_ses-*_bold.json',
                  'description': 'blood-oxygen-level-dependent (BOLD) JSON file'}

""" PET """

# BIDS

PET_FDG_NII = {'pattern': 'pet/sub-*_ses-*_acq-fdg_pet.nii*',
               'description': 'FDG-PET data'}

PET_FDG_JSON = {'pattern': 'pet/sub-*_ses-*_acq-fdg_pet.json',
                'description': 'JSON file describing the point spread function (PSF) in FDG PET.'}

PET_AV45_NII = {'pattern': 'pet/sub-*_ses-*_acq-av45_pet.nii*',
                'description': 'AV45-PET data'}

PET_AV45_JSON = {'pattern': 'pet/sub-*_ses-*_acq-av45_pet.json',
                 'description': 'JSON file describing the point spread function (PSF) in AV45 PET.'}
