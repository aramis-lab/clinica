# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 15:56:18 2016

@author: jacquemont
"""
from clinica.pipeline.preproc.DWI_launch_preproc import launch
import os
from os.path import realpath,split,join

data_path = join(split(realpath(__file__))[0], 'data/DWI_launch_preproc')

DWI = join(data_path, 'DWI.nii')
T1 = join(data_path, 'T1.nii')
b_values = join(data_path, 'b_values.txt')
b_vectors = join(data_path, 'b_vectors.txt')

working_direct = os.getcwd()
datasink_direct = os.getcwd()

launch(DWI, T1, b_values, b_vectors, working_direct, datasink_direct)