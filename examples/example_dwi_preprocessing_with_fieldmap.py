#!/usr/bin/python

from __future__ import absolute_import

from clinica.pipeline.dwi.dwi_preprocessing import diffusion_preprocessing_fieldmap_based
from clinica.pipeline.dwi.dwi_preprocessing_utils import count_b0s
import nipype.interfaces.fsl as fsl

import os
from os.path import realpath,split,join
import tempfile


data_path = join(split(realpath(__file__))[0], 'external-data/BIDS-example/sub-CLNC01/ses-M00/')
output_directory = tempfile.mkdtemp()

print("Datasink Directory -> %s" % output_directory)

preprocessing = diffusion_preprocessing_fieldmap_based(subject_id='sub-CLNC01', session_id='ses-M00',
                                                       analysis_series_id='default', caps_directory=output_directory,
                                                       delta_te=2.46e-3, echo_spacing=0.39e-3,
                                                       num_b0s=count_b0s(join(data_path, 'dwi/sub-CLNC01_ses-M00_dwi.bval'))
                                                       )
preprocessing.inputs.inputnode.in_dwi   = join(data_path, 'dwi/sub-CLNC01_ses-M00_dwi.nii.gz')
preprocessing.inputs.inputnode.in_bvals  = join(data_path, 'dwi/sub-CLNC01_ses-M00_dwi.bval')
preprocessing.inputs.inputnode.in_bvecs  = join(data_path, 'dwi/sub-CLNC01_ses-M00_dwi.bvec')
preprocessing.inputs.inputnode.in_fmap_mag  = join(data_path, 'fmap/sub-CLNC01_ses-M00_magnitude1.nii.gz')
preprocessing.inputs.inputnode.in_fmap_pha  = join(data_path, 'fmap/sub-CLNC01_ses-M00_phasediff.nii.gz')

preprocessing.run()

print("Datasink Directory -> %s" % output_directory)
