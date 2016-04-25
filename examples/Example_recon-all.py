#!/usr/bin/python#
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 15:56:18 2016

@author: Junhao WEN
"""

from __future__ import absolute_import
from clinica.pipeline.preprocessing.t1_freesurfer_workflow import recon_all_pipeline
import os
from os.path import realpath,split,join
import tempfile
from nipype.interfaces.freesurfer.preprocess import ReconAll


try:
    if ReconAll.version.fget.func_globals['__version__'].split(".") < ['0','11','0']:
        raise RuntimeError('ReconAll version should at least be version of 0.11.0')
except Exception as e:
    print(str(e))
    exit(1)

experiment_dir = join(split(realpath(__file__))[0], 'data')

data_dir = join(experiment_dir, 'Recon-all')
output_dir = tempfile.mkdtemp()
working_dir = tempfile.mkdtemp()

print("Data Directory -> %s" % data_dir)
print("Output Directory -> %s" % output_dir)

print("Running...")
T1_recon_all= recon_all_pipeline(data_dir, experiment_dir, output_dir, working_dir)
T1_recon_all.run("MultiProc", plugin_args={'n_procs':4})

