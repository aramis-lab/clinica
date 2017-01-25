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

data_path = join(split(realpath(__file__))[0], 'external-data/DWI_postproc')

in_parcellation = join('path to parc')
configuration_file = join(data_path, 'config_file.nii')
atlas_scalar_image = join(data_path, 'WM_atlas_FA.nii')
lut_type = 'lut_type'
lut_path = join('path to LUT')
in_tracks = join(data_path, 'Tracks.tck')
connectome_metric = 'count'
zeros_diagonal=True

working_directory = tempfile.mkdtemp()
datasink_directory = tempfile.mkdtemp()

print("Working Directory -> %s" % working_directory)
print("Datasink Directory -> %s" % datasink_directory)

print("Running connectome construction")
connectome =  connectome_construction_pipeline(in_parcellation, configuration_file, lut_type, lut_path, in_tracks, connectome_metric, working_directory, datasink_directory, in_scalar_image='', zeros_diagonal=True)
connectome.run()
print("connectome construction done")

print("Working Directory -> %s" % working_directory)
print("Datasink Directory -> %s" % datasink_directory)

