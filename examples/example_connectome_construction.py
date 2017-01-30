# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 15:56:18 2016

@author: jacquemont
"""
from __future__ import absolute_import
from clinica.pipeline.dwi.dwi_connectome_construction import connectome_construction_pipeline
import os
from os.path import realpath,split,join
import tempfile

data_path = join(split(realpath(__file__))[0], 'external-data/CLNC01_M00')

in_parcellation = join(data_path, 'mri', 'aparc.a2009s+aseg.mgz')
configuration_file = join(split(realpath(__file__))[0], 'config_fs_a2009s.txt')
lut_type = 'freesurfer'
lut_path = join('/Applications/freesurfer/FreeSurferColorLUT.txt')
in_tracks = join('/tmp/subject_example/analysis-series-default/subjects/CLNC01/M00/dwi/mrtrix/sub-CLNC01_ses-M00_fibers-100K.tck')
connectome_metric = 'count'
zero_diagonal=True

working_directory = tempfile.mkdtemp()
datasink_directory = tempfile.mkdtemp()

print("Working Directory -> %s" % working_directory)
print("Datasink Directory -> %s" % datasink_directory)

print("Running connectome construction")
connectome =  connectome_construction_pipeline(
    in_parcellation=in_parcellation, configuration_file=configuration_file, lut_type=lut_type, lut_path=lut_path, in_tracks=in_tracks, connectome_metric=connectome_metric,
    working_directory=working_directory, datasink_directory=datasink_directory,in_scalar_image=None, zero_diagonal=zero_diagonal
)

#connectome_construction_pipeline(in_parcellation, configuration_file, lut_type, lut_path, in_tracks, connectome_metric, working_directory, datasink_directory, in_scalar_image='', zeros_diagonal=True)
connectome.run()
print("connectome construction done")

print("Working Directory -> %s" % working_directory)
print("Datasink Directory -> %s" % datasink_directory)

