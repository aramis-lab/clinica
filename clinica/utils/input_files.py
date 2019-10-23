# coding: utf8

"""
Describe files to grab, to use with inputs.clinica_file_reader() and inputs.clinica_group_reader()

"""

""" T1w """

# BIDS

T1W_NII = {'pattern': 'sub-*_ses-*_t1w.nii*',
           'description': 'T1w MRI'}

# FreeSurfer

T1_FS_WM = {'pattern': 't1/freesurfer_cross_sectional/sub-*_ses-*/mri/wm.seg.mgz',
            'description': 'segmentation of white matter (mri/wm.seg.mgz)',
            'needed_pipeline': 't1-freesurfer'}

T1_FS_BRAIN = {'pattern': 't1/freesurfer_cross_sectional/sub-*_ses-*/mri/brain.mgz',
               'description': ' extracted brain from T1w MRI',
               'needed_pipeline': 't1-freesurfer'}

T1_FS_ORIG_NU = {'pattern': 't1/freesurfer_cross_sectional/sub-*_ses-*/mri/orig_nu.mgz',
                 'description': 'intensity normalized volume generated after correction for'
                                ' non-uniformity in FreeSurfer (file mri/orig_nu.mgz) ',
                 'needed_pipeline': 't1-freesurfer'}

T1_FS_WM_SURF_R = {'pattern': 't1/freesurfer_cross_sectional/sub-*_ses-*/surf/rh.white',
                   'description': 'right white matter/gray matter border surface.',
                   'needed_pipeline': 't1-freesurfer'}

T1_FS_WM_SURF_L = {'pattern': 't1/freesurfer_cross_sectional/sub-*_ses-*/surf/lh.white',
                   'description': 'left white matter/gray matter border surface.',
                   'needed_pipeline': 't1-freesurfer'}

T1_FS_DESTRIEUX = {'pattern': 't1/freesurfer_cross_sectional/sub-*_ses-*/mri/aparc.a2009s+aseg.mgz',
                   'description': 'Destrieux-based segmentation (mri/aparc.a2009s+aseg.mgz)',
                   'needed_pipeline': 't1-freesurfer'}

T1_FS_DESTRIEUX_PARC_L = {'pattern': 't1/freesurfer_cross_sectional/sub-*_ses-*/label/lh.aparc.a2009s.annot',
                          'description': 'left hemisphere surface-based Destrieux parcellation (label/lh.aparc.a2009s.annot)',
                          'needed_pipeline': 't1-freesurfer'}

T1_FS_DESTRIEUX_PARC_R = {'pattern': 't1/freesurfer_cross_sectional/sub-*_ses-*/label/rh.aparc.a2009s.annot',
                          'description': 'right hemisphere surface-based Destrieux parcellation (label/rh.aparc.a2009s.annot)',
                          'needed_pipeline': 't1-freesurfer'}

T1_FS_DESIKAN = {'pattern': 't1/freesurfer_cross_sectional/sub-*_ses-*/mri/aparc+aseg.mgz',
                 'description': 'Desikan-based segmentation (mri/aparc.a2009s+aseg.mgz)',
                 'needed_pipeline': 't1-freesurfer'}

T1_FS_DESIKAN_PARC_L = {'pattern': 't1/freesurfer_cross_sectional/sub-*_ses-*/label/lh.aparc.annot',
                        'description': 'left hemisphere surface-based Desikan parcellation (label/lh.aparc.annot)',
                        'needed_pipeline': 't1-freesurfer'}

T1_FS_DESIKAN_PARC_R = {'pattern': 't1/freesurfer_cross_sectional/sub-*_ses-*/label/rh.aparc.annot',
                        'description': 'right hemisphere surface-based Desikan parcellation (label/rh.aparc.annot)',
                        'needed_pipeline': 't1-freesurfer'}

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
                         'description': 'b0 brainmaisk',
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
