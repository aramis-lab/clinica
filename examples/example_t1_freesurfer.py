#!/usr/bin/python#
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 15:04:17 2016

@author: Junhao WEN

Example to lanch freesurfer-recon-all pipeline
"""

from __future__ import absolute_import
from clinica.pipeline.t1.t1_freesurfer_workflows import t1_freesurfer_pipeline
import os
from os.path import realpath,split,join
import tempfile, errno

output_dir = tempfile.mkdtemp()
data_path = join(split(realpath(__file__))[0], 'external-data/BIDS-example/sub-CLNC01/ses-M00/')
anat_t1 = os.path.join(data_path, 'anat/sub-CLNC01_ses-M00_T1w.nii.gz')
subject_list = 'sub-CLNC01'
session_list = 'ses-M00'
subjects_dir = os.path.join(output_dir + 'analysis-series-default/subjects/' + subject_list + '/' + session_list + '/' + 't1' + '/' + 'freesurfer-cross-sectional')
try:
    os.makedirs(subjects_dir)
except OSError as exception:
    if exception.errno != errno.EEXIST:
        raise


freesurfer_t1 = t1_freesurfer_pipeline(output_dir, 'default', output_dir, '-qcache')

freesurfer_t1.inputs.recon_all.subjects_dir = subjects_dir
freesurfer_t1.inputs.recon_all.subject_id = subject_list + '_' + session_list
freesurfer_t1.inputs.recon_all.T1_files = anat_t1
freesurfer_t1.inputs.flagnode.t1_list = anat_t1
freesurfer_t1.inputs.lognode.subject_list = subject_list
freesurfer_t1.inputs.lognode.session_list = session_list

print("Results will be stored in the following path: %s" % output_dir)
freesurfer_t1.run()
print("Results are stored here: %s" % output_dir)

