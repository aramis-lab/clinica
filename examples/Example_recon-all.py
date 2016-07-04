#!/usr/bin/python#
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 15:04:17 2016

@author: Junhao WEN
"""

from __future__ import absolute_import
from clinica.pipeline.preprocessing.t1_freesurfer_workflow import recon_all_pipeline
import os.path
import tempfile

# data_path = join(split(realpath(__file__))[0], 'data/DWI_postproc')

data_dir = os.path.abspath('data/Recon-all')
output_dir = tempfile.mkdtemp()

print("Data Directory -> %s" % data_dir)
print("Output Directory -> %s" % output_dir)
print("Running...")


def example_with_structure_A():
    return recon_all_pipeline('%s/%s/x.nii', output_dir, 3)

def example_with_structure_B():
    return recon_all_pipeline('%s/struct/x.nii', output_dir, 3)

T1_recon_all= recon_all_pipeline(data_dir, output_dir, 3)
# T1_recon_all= recon_all_pipeline(['data/%s/f%d.nii', 'extra/d1/m1.nii', 'more/x.nii'], output_dir, 3)
# T1_recon_all= recon_all_pipeline('%s/struct.nii', output_dir, 3)
# T1_recon_all= recon_all_pipeline(['data/example1/m1.nii', 'data/example2/file.nii', 'more/f1.nii', 'nore/f2.nii'], output_dir, 3)
T1_recon_all.run("MultiProc", plugin_args={'n_procs':4})
