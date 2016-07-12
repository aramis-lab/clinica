# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 15:56:18 2016

@author: jacquemont
"""
from __future__ import absolute_import
from clinica.pipeline.preprocessing.dwi_launch_preproc import create_dwi_preproc_syb
from clinica.pipeline.preprocessing.dwi_utils import count_b0s
import os
from os.path import realpath,split,join
import tempfile

import nipype.interfaces.fsl as fsl

try:
    if fsl.Info.version().split(".") < ['5','0','5']:
        raise RuntimeError('FSL version must be great then 5.0.5')
except Exception as e:
    print(str(e))
    exit(1)

data_path = join(split(realpath(__file__))[0], 'external-data/DWI_launch_preproc')

DWI = join(data_path, 'DWI.nii')
T1 = join(data_path, 'T1.nii')
b_values = join(data_path, 'b_values.txt')
b_vectors = join(data_path, 'b_vectors.txt')

working_direct = tempfile.mkdtemp()
datasink_direct = tempfile.mkdtemp()

print("Working Directory -> %s" % working_direct)
print("Datasink Directory -> %s" % datasink_direct)

print("Running...")
dwi_preproc_syb = create_dwi_preproc_syb(DWI, T1, b_values, b_vectors, working_direct, count_b0s(b_values), datasink_direct)
dwi_preproc_syb.run()

print("Working Directory -> %s" % working_direct)
print("Datasink Directory -> %s" % datasink_direct)
