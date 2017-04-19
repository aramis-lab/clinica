#!/usr/bin/python

from __future__ import absolute_import

from clinica.pipeline.dwi.dwi_preprocessing import diffusion_preprocessing_phasediff_fieldmap
from clinica.pipeline.dwi.dwi_preprocessing_utils import count_b0s
from os.path import realpath,split,join
import tempfile


path_to_data = join(split(realpath(__file__))[0], 'external-data/BIDS-example/sub-CLNC01/ses-M00/')
output_directory = tempfile.mkdtemp()

print("Datasink Directory -> %s" % output_directory)

preprocessing = diffusion_preprocessing_phasediff_fieldmap(
    participant_id='sub-CLNC01', session_id='ses-M00',
    caps_directory=output_directory,
    delta_echo_time=2.46e-3,
    effective_echo_spacing=0.39e-3,
    phase_encoding_direction='y',
    num_b0s=count_b0s(join(path_to_data, 'dwi/sub-CLNC01_ses-M00_dwi.bval'))
    )
preprocessing.inputs.inputnode.in_dwi             = join(path_to_data, 'dwi/sub-CLNC01_ses-M00_dwi.nii.gz')
preprocessing.inputs.inputnode.in_bvals           = join(path_to_data, 'dwi/sub-CLNC01_ses-M00_dwi.bval')
preprocessing.inputs.inputnode.in_bvecs           = join(path_to_data, 'dwi/sub-CLNC01_ses-M00_dwi.bvec')
preprocessing.inputs.inputnode.in_fmap_magnitude  = join(path_to_data, 'fmap/sub-CLNC01_ses-M00_magnitude1.nii.gz')
preprocessing.inputs.inputnode.in_fmap_phasediff  = join(path_to_data, 'fmap/sub-CLNC01_ses-M00_phasediff.nii.gz')

preprocessing.run()

print("Datasink Directory -> %s" % output_directory)
