#!/usr/bin/python

"""This module launches the T1-FSL pipeline."""

from __future__ import absolute_import
from clinica.pipeline.t1.t1_fsl import t1_fsl_segmentation_pipeline

from os.path import realpath,split,join
import tempfile

path_to_data = join(split(realpath(__file__))[0], 'external-data/BIDS-example/sub-CLNC01/ses-M00/')
output_directory = tempfile.mkdtemp()


fsl_t1 = t1_fsl_segmentation_pipeline(
    participant_id='sub-CLNC01', session_id='ses-M00',
    caps_directory=output_directory, working_directory=None,
    is_bias_corrected=False)

fsl_t1.inputs.inputnode.in_t1 = join(path_to_data, 'anat/sub-CLNC01_ses-M00_T1w.nii.gz')

print("Results will be stored in the following path: %s" % output_directory)
fsl_t1.run()
print("Results are stored here: %s" % output_directory)
