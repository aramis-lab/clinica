# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 17:03:04 2016

@author: wen
"""
from __future__ import absolute_import
from clinica.pipeline.preprocessing.SurfStatPipeline import SurfStatGruopAnalysis
from os.path import realpath,split,join

a_required_path = join(split(realpath(__file__))[0], 'ClinicaSurfstat');
PATH_TO_RECON_ALL_OUTPUTS = join(split(realpath(__file__))[0], 'data/Recon-all_Output');
CSVFilename  = join(a_required_path, 'Database/template.csv')
Format = '%s %s %s %f';
ContrastLinearModel = '1 + Label + Gender + Age';
SurfStatGruopAnalysis(ContrastLinearModel,Format, CSVFilename, PATH_TO_RECON_ALL_OUTPUTS, a_required_path)