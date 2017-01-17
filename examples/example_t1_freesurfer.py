#!/usr/bin/python#
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 15:04:17 2016

@author: Junhao WEN
"""

from __future__ import absolute_import
from clinica.pipeline.t1.t1_freesurfer import datagrabber_t1_freesurfer_pipeline
from os.path import realpath,split,join
import time
import tempfile

BIDS_dir ='/Volumes/dataARAMIS/users/CLINICA/CLINICA_datasets/BIDS/PREVDEMALS_BIDS/GENFI'
output_dir ='~/test/test-reconall'
tsv_file ='/Volumes/dataARAMIS/users/junhao.wen/PhD/PREVDEMALS/Freesurfer/Reconall/reconall_GENFI/clinica_reconall_result/prevdemals_67subjs/subjects_sessions_list_dwitest.tsv'
start = time.time()
# working_directory='~/test/test-reconall'

# this is the example to run CAPP dataset
def recon_all_example_CAPP():
    return datagrabber_t1_freesurfer_pipeline(BIDS_dir, output_dir, tsv_file)

if 2 > 1:
    print("Data Directory -> %s" % BIDS_dir)
    print("Output Directory -> %s" % output_dir)
    print("Running...")
    T1_recon_all = recon_all_example_CAPP()
    T1_recon_all.run("MultiProc", plugin_args={'n_procs':4})
    time_consuming = time.time() - start
    print 'END! time consuming is : %s' % time_consuming
#
# from clinica.pipeline.t1.t1_freesurfer_utils import get_dirs_check_reconalled
#
# get_dirs_check_reconalled('~/test/test-reconall', '/Volumes/dataARAMIS/users/junhao.wen/PhD/PREVDEMALS/Freesurfer/Reconall/reconall_GENFI/clinica_reconall_result/prevdemals_67subjs/subjects_sessions_list_dwitest.tsv', 'default' )



    # command line example:
# clinica run t1-freesurfer /aramis/dataARAMIS/users/CLINICA/CLINICA_datasets/BIDS/PREVDEMALS_BIDS/GENFI ~/test/test-reconall-lab/ /aramis/dataARAMIS/users/junhao.wen/PhD/PREVDEMALS/Freesurfer/Reconall/reconall_GENFI/clinica_reconall_result/subjects_visits_list_PREVDEMALS.tsv 'default'