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
temporary_dir = tempfile.mkdtemp()
output_dir = tempfile.mkdtemp()


print("Data Directory -> %s" % data_dir)
print("Output Directory -> %s" % output_dir)
print("temporary ouput Directory -> %s" % temporary_dir)

print("Running...")
T1_recon_all= recon_all_pipeline(data_dir,temporary_dir, output_dir)
T1_recon_all.run("MultiProc", plugin_args={'n_procs':4})

