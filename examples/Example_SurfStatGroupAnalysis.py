# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 17:03:04 2016

@author: wen
"""
from __future__ import absolute_import
from clinica.pipeline.preprocessing.SurfStatPipeline import SurfStatGruopAnalysis
#from os.path import realpath,split,join
import os

#a_required_path = join(split(realpath(__file__))[0], 'ClinicaSurfstat');
path = os.getcwd()
parent_path = os.path.dirname(path)
a_required_path = os.path.join(parent_path, 'lib/ClinicaSurfstat')
PATH_TO_RECON_ALL_OUTPUTS = os.path.join(os.path.split(os.path.realpath(__file__))[0], 'data/Recon-all_Output');
CSVFilename  = os.path.join(os.path.split(os.path.realpath(__file__))[0], 'ClinicaSurfstat/Database/template.csv');
Format = '%s %s %s %f';
ContrastLinearModel = '1 + Label + Gender + Age';
result_dir = os.path.join(os.path.split(os.path.realpath(__file__))[0], 'ClinicaSurfstat');
SurfStat = SurfStatGruopAnalysis(ContrastLinearModel,Format, CSVFilename, PATH_TO_RECON_ALL_OUTPUTS, a_required_path, result_dir)
SurfStat.run()