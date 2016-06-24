#!/usr/bin/python#
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 15:04:17 2016

@author: Junhao WEN
"""

from __future__ import absolute_import
from clinica.pipeline.preprocessing.t1_freesurfer_workflow import recon_all_pipeline
import os.path
import tmpfile

data_dir = os.path.abspath('data/Recon-all')
output_dir = tmpfile.mkdtemp()

print("Data Directory -> %s" % data_dir)
print("Output Directory -> %s" % output_dir)
print("Running...")

T1_recon_all= recon_all_pipeline(data_dir, output_dir, 3)
T1_recon_all.run("MultiProc", plugin_args={'n_procs':4})
