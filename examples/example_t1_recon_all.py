#!/usr/bin/python#
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 15:04:17 2016

@author: Junhao WEN
"""

from __future__ import absolute_import
from clinica.pipeline.t1.t1_freesurfer import recon_all_pipeline
from os.path import realpath,split,join
import tempfile

# Test for CAPP dataset
data_dir = join(split(realpath(__file__))[0], 'external-data/CAPP_BIDSsource')
output_dir = tempfile.mkdtemp()
tsv_file = join(split(realpath(__file__))[0], 'external-data/subjects_visits_list_CAPP.tsv')
dataset_name = 'CAPP'

print("Data Directory -> %s" % data_dir)
print("Output Directory -> %s" % output_dir)
print("Running...")

# this is the example to run CAPP dataset
def recon_all_example_CAPP():
    return recon_all_pipeline(data_dir, output_dir,tsv_file, dataset_name)

# this is the example to run INSIGHT dataset
def recon_all_example_INSIGHT():
    return recon_all_pipeline(data_dir, output_dir,tsv_file, dataset_name)

if dataset_name == 'CAPP':
    T1_recon_all = recon_all_example_CAPP()
    T1_recon_all.run("MultiProc", plugin_args={'n_procs':4})
else:
    T1_recon_all = recon_all_example_INSIGHT()
    T1_recon_all.run("MultiProc", plugin_args={'n_procs':4})

