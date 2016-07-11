# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 17:03:04 2016
    Explaination: This example is to show you how to use ClinicaSurfStat to do group analysis, basically, the parameters below:
      :param: ContrastLinearModel: string, the linear model that fit into the GLM, for example '1+Lable'.
      :param: Format: string, the format that you want to use for your CSV file column variables, it depends on the CSV file. 
      :param: CSVFilename: string, the path to your csv file, we put it in the directory 'ClinicaSurfstat/Database'
      :param: PATH_TO_RECON_ALL_OUTPUTS:  the output file from recon-all pipeline,specifically, files: ?h.thickness.fwhm**.mgh. 
              we put in the directory 'data/Recon-all_Output'.
      :param: a_required_path:  this is the path to find the matlabscript(NipypeSurfStat.m) that you are going to run.
      :param: result_dir: the directory to contain the result images. The out put directory in the example is in the directory 'ClinicaSurfstat/Figures   '.

@author: wen
"""
from __future__ import absolute_import
from clinica.pipeline.preprocessing.SurfStatPipeline import SurfStatGruopAnalysis
from os.path import realpath,split,join
import os

path = join(split(realpath(__file__))[0])
parent_path = os.path.dirname(path)
a_required_path = os.path.join(parent_path, 'lib/ClinicaSurfstat')
PATH_TO_RECON_ALL_OUTPUTS = os.path.join(os.path.split(os.path.realpath(__file__))[0], 'data/Recon-all_Output')
CSVFilename  = os.path.join(os.path.split(os.path.realpath(__file__))[0], 'ClinicaSurfstat/Database/template.csv')
Format = '%s %s %s %f'
ContrastLinearModel = '1 + Label + Gender + Age'
result_dir = os.path.join(os.path.split(os.path.realpath(__file__))[0], 'ClinicaSurfstat')
SurfStat = SurfStatGruopAnalysis(ContrastLinearModel,Format, CSVFilename, PATH_TO_RECON_ALL_OUTPUTS, a_required_path, result_dir)
SurfStat.run()