#!/usr/bin/python#
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 15:04:17 2016

@author: Junhao WEN
"""

from __future__ import absolute_import
from clinica.pipeline.t1.t1_freesurfer import t1_freesurfer_tsv_pipeline
from os.path import realpath,split,join
import time
import tempfile

output_dir = '/Volumes/dataARAMIS/users/junhao.wen/PhD/PREVDEMALS/Freesurfer/Reconall/reconall_GENFI/clinica_reconall_result/micheal_computer'
tsv_file_CAPP = '/Volumes/dataARAMIS/users/junhao.wen/PhD/PREVDEMALS/Freesurfer/Reconall/reconall_GENFI/clinica_reconall_result/subjects_visits_list_PREVDEMALS.tsv'

start = time.time()
working_directory='~/test/test-reconall-lab'
analysis_series_id = 'michael'
T1_recon_all_statistics = t1_freesurfer_tsv_pipeline(output_dir, tsv_file_CAPP, analysis_series_id = analysis_series_id, working_directory=working_directory)
T1_recon_all_statistics.run("MultiProc", plugin_args={'n_procs':4})
time_consuming = time.time() - start
print 'END! time consuming is : %s' % time_consuming