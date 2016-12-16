# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 15:56:18 2016

@author: jacquemont
"""
from clinica.pipeline.dwi.dwi_preprocessing import create_dwi_preproc_syb
from clinica.pipeline.dwi.dwi_preprocessing_utils import count_b0s

import os
from os.path import realpath,split,join
import tempfile


data_path = join(split(realpath(__file__))[0], 'external-data/DWI_launch_preproc')

DWI = join(data_path, 'DWI.nii')
T1 = join(data_path, 'T1.nii')
b_values = join(data_path, 'b_values.txt')
b_vectors = join(data_path, 'b_vectors.txt')

number_of_b0s = count_b0s(b_values)
print("Number of B0 -> %d" % number_of_b0s)


working_direct = tempfile.mkdtemp()
datasink_direct = tempfile.mkdtemp()
print("Working Directory -> %s" % working_direct)
print("Datasink Directory -> %s" % datasink_direct)

print("Running...")
dwi_preproc_syb = create_dwi_preproc_syb(DWI, T1, b_values, b_vectors, working_direct, number_of_b0s, datasink_direct)
dwi_preproc_syb.run('MultiProc', plugin_args={'n_procs': 4})

print("Working Directory -> %s" % working_direct)
print("Datasink Directory -> %s" % datasink_direct)
