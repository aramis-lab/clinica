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

data_path = join(split(realpath(__file__))[0], 'external-data/raw_data/subject_example')

output_directory = tempfile.mkdtemp()

print("Datasink Directory -> %s" % output_directory)

preprocessing = diffusion_preprocessing_fieldmap_based(datasink_directory=output_directory,
    num_b0s=count_b0s(join(data_path, 'DWI/subject_example_dwi.bval')))

preprocessing.inputs.inputnode.in_file   = join(data_path, 'DWI/subject_example_dwi.nii.gz')
preprocessing.inputs.inputnode.in_bvals  = join(data_path, 'DWI/subject_example_dwi.bval')
preprocessing.inputs.inputnode.in_bvecs  = join(data_path, 'DWI/subject_example_dwi.bvec')
preprocessing.inputs.inputnode.bmap_mag  = join(data_path, 'B0MAP/subject_example_b0map_echo1.nii.gz')
preprocessing.inputs.inputnode.bmap_pha  = join(data_path, 'B0MAPph/subject_example_b0mapph.nii.gz')

preprocessing.run()

print("Datasink Directory -> %s" % output_directory)
