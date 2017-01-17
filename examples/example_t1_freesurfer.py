#!/usr/bin/python#
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 15:04:17 2016

@author: Junhao WEN

Example to lanch freesurfer-recon-all pipeline
"""

from __future__ import absolute_import
from clinica.pipeline.t1.t1_freesurfer import datagrabber_t1_freesurfer_pipeline
import time

BIDS_dir ='path to your BIDS dataset'
output_dir ='path to your CAPS dataset'
start = time.time()

# this is the example to run CAPP dataset
def recon_all_example_CAPP():
    return datagrabber_t1_freesurfer_pipeline(BIDS_dir, output_dir)

if 2 > 1:
    print("Data Directory -> %s" % BIDS_dir)
    print("Output Directory -> %s" % output_dir)
    print("Running...")
    T1_recon_all = recon_all_example_CAPP()
    T1_recon_all.run("MultiProc", plugin_args={'n_procs':4})
    time_consuming = time.time() - start
    print 'END! time consuming is : %s' % time_consuming

    # command line example:
# clinica run t1-freesurfer /aramis/dataARAMIS/users/CLINICA/CLINICA_datasets/BIDS/PREVDEMALS_BIDS/GENFI ~/test/test-reconall-lab/ /aramis/dataARAMIS/users/junhao.wen/PhD/PREVDEMALS/Freesurfer/Reconall/reconall_GENFI/clinica_reconall_result/subjects_visits_list_PREVDEMALS.tsv 'default'