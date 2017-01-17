#!/usr/bin/python

from __future__ import absolute_import

from clinica.pipeline.dwi.dwi_preprocessing import diffusion_preprocessing_fieldmap_based
from clinica.pipeline.dwi.dwi_preprocessing_utils import count_b0s
import nipype.interfaces.fsl as fsl

import os
from os.path import realpath,split,join
import tempfile

try:
    if fsl.Info.version() is None or fsl.Info.version().split(".") < ['5','0','5']:
        raise RuntimeError('FSL version must be great then 5.0.5')
except Exception as e:
    print(str(e))
    exit(1)

data_path = join(split(realpath(__file__))[0], 'external-data/BIDS-example/sub-CLNC01/ses-M00/')

output_directory = tempfile.mkdtemp()

print("Datasink Directory -> %s" % output_directory)

preprocessing = diffusion_preprocessing_fieldmap_based(subject_id='sub-CLNC01', session_id='ses-M00',
                                                       analysis_series_id='default', caps_directory=output_directory,
                                                       num_b0s=count_b0s(join(data_path, 'dwi/sub-CLNC01_ses-M00_dwi.bval'))
                                                       )
preprocessing.inputs.inputnode.in_file   = join(data_path, 'dwi/sub-CLNC01_ses-M00_dwi.nii.gz')
preprocessing.inputs.inputnode.in_bvals  = join(data_path, 'dwi/sub-CLNC01_ses-M00_dwi.bval')
preprocessing.inputs.inputnode.in_bvecs  = join(data_path, 'dwi/sub-CLNC01_ses-M00_dwi.bvec')
preprocessing.inputs.inputnode.bmap_mag  = join(data_path, 'fmap/sub-CLNC01_ses-M00_magnitude1.nii.gz')
preprocessing.inputs.inputnode.bmap_pha  = join(data_path, 'fmap/sub-CLNC01_ses-M00_phasediff.nii.gz')

preprocessing.run()

print("Datasink Directory -> %s" % output_directory)
