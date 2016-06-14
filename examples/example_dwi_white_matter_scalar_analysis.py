# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 15:56:18 2016

@author: jacquemont
"""
from __future__ import absolute_import
from clinica.pipeline.postprocessing.dwi_white_matter_scalar_analysis import create_connection_wise_analysis_pipeline
import os
from os.path import realpath,split,join
import tempfile

data_path = join(split(realpath(__file__))[0], 'data/DWI_postproc')

in_scalar_image = join(data_path, 'FA.nii')
atlas_labels = join(data_path, 'WM_atlas_labels.nii')
atlas_scalar_image = join(data_path, 'WM_atlas_FA.nii')

working_direct = tempfile.mkdtemp()
datasink_direct = tempfile.mkdtemp()

print("Working Directory -> %s" % working_direct)
print("Datasink Directory -> %s" % datasink_direct)

print("Running DTI atlas scalar analysis")
DTI_atlas_scalar_analysis =  create_DTI_atlas_scalar_analysis(in_scalar_image, atlas_labels, atlas_scalar_image, working_directory, datasink_directory)
DTI_atlas_scalar_analysis.run()
print("DTI atlas scalar analysis done")

print("Working Directory -> %s" % working_direct)
print("Datasink Directory -> %s" % datasink_direct)

