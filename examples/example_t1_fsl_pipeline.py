#!/usr/bin/python

"""This module launches the T1-FSL pipeline."""

from __future__ import absolute_import
from clinica.pipeline.t1.t1_fsl import t1_fsl_segmentation_pipeline


from os.path import realpath,split,join
import tempfile


data_path = join(split(realpath(__file__))[0], 'external-data/raw_data/subject_example')
output_directory = tempfile.mkdtemp()


fsl_t1 = t1_fsl_segmentation_pipeline(caps_directory=output_directory, is_bias_corrected=False)
fsl_t1.inputs.inputnode.in_t1 = join(data_path, 'T1/subject_example_t1.nii.gz')


print("Datasink Directory -> %s" % output_directory)
fsl_t1.run()
print("Datasink Directory -> %s" % output_directory)
