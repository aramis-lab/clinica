#!/usr/bin/python

from __future__ import absolute_import

from clinica.pipeline.preprocessing.dwi_launch_preproc import diffusion_preprocessing_fieldmap_based
import nipype.interfaces.fsl as fsl

import os
from os.path import realpath,split,join
import tempfile

try:
    if fsl.Info.version().split(".") < ['5','0','5']:
        raise RuntimeError('FSL version must be great then 5.0.5')
except Exception as e:
    print(str(e))
    exit(1)

data_path = join(split(realpath(__file__))[0], 'data/raw_data/subject_example')

output_directory = "/tmp/output_dwi_preprocessing_fieldmap_based/new_suject"

preprocessing = diffusion_preprocessing_fieldmap_based(datasink_directory=output_directory)

preprocessing.inputs.inputnode.in_file   = join(data_path, 'DWI/subject_example_dwi.nii.gz')
preprocessing.inputs.inputnode.in_bvals  = join(data_path, 'DWI/subject_example_dwi.bval')
preprocessing.inputs.inputnode.in_bvecs  = join(data_path, 'DWI/subject_example_dwi.bvec')
preprocessing.inputs.inputnode.bmap_mag  = join(data_path, 'B0MAP/0000.nii.gz')
preprocessing.inputs.inputnode.bmap_pha  = join(data_path, 'B0MAPph/subject_example_b0mapph.nii.gz')



#working_direct = tempfile.mkdtemp()
#datasink_direct = tempfile.mkdtemp()

#print("Working Directory -> %s" % working_direct)
#print("Datasink Directory -> %s" % datasink_direct)

print("Running...")
preprocessing.run()
