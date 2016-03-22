# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 15:56:18 2016

@author: jacquemont
"""
from clinica.pipeline.preproc.DWI_launch_preproc import launch
import os
from os.path import abspath,join

DWI = join(abspath(__file__), '/data/DWI_launch_preproc/DWI.nii')
T1 = join(abspath(__file__), '/data/DWI_launch_preproc/T1.nii')
b_values = join(abspath(__file__), '/data/DWI_launch_preproc/b_values.txt')
b_vectors = join(abspath(__file__), '/data/DWI_launch_preproc/b_vectors.txt')

working_direct = os.getcwd()
datasink_direct = os.getcwd()

launch(DWI, T1, b_values, b_vectors, working_direct, datasink_direct)