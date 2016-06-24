#!/usr/bin/python#
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 15:04:17 2016

@author: Junhao WEN
"""

from __future__ import absolute_import
from clinica.pipeline.preprocessing.t1_freesurfer_workflow import recon_all_pipeline
from os.path import realpath,split,join
import tempfile

data_dir = join(split(realpath(__file__))[0], 'data/Recon-all')
field_templateA = '%s/%s.nii'
template_argsA =  [['subject_id', 'struct']]
field_templateB = '%s/%s/*.nii'
template_argsB = [['subject_id', 'dir']]
output_dir = tempfile.mkdtemp()

print("Data Directory -> %s" % data_dir)
print("Output Directory -> %s" % output_dir)
print("Running...")

def recon_all_exampleA():
    return recon_all_pipeline(data_dir, output_dir, 3, field_templateA, template_argsA)
    
def recon_all_exampleB():
    return recon_all_pipeline(data_dir, output_dir, 3, field_templateB, template_argsB)
T1_recon_allA = recon_all_exampleA()
T1_recon_allA.run("MultiProc", plugin_args={'n_procs':4})
T1_recon_allB = recon_all_exampleB()
T1_recon_allB.run("MultiProc", plugin_args={'n_procs':4})
