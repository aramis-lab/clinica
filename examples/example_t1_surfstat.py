# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 17:03:04 2016
    Explaination: This example is to show you how to use Clinicasurfstat to do group analysis, basically, the parameters below:
      :param: ContrastLinearModel: string, the linear model that fit into the GLM, for example '1+Lable'.
      :param: Format: string, the format that you want to use for your CSV file column variables, it depends on the CSV file.
      :param: CSVFilename: string, the path to your csv file, we put it in the directory 'Clinicasurfstat/Database'
      :param: PATH_TO_RECON_ALL_OUTPUTS:  the output file from recon-all pipeline,specifically, files: ?h.thickness.fwhm**.mgh.
              we put in the directory 'data/Recon-all_Output'.
      :param: a_required_path:  this is the path to find the matlabscript(Nipypesurfstat.m) that you are going to run.
      :param: result_dir: the directory to contain the result images. The out put directory in the example is in the directory 'Clinicasurfstat/Figures   '.

@author: wen
"""
from __future__ import absolute_import
from clinica.pipeline.postprocessing.t1_surfstat_workflow import clinica_surfstat
from os.path import realpath, split, join
import tempfile
import time


#path = join(split(realpath(__file__))[0])
#parent_path = dirname(dirname(path))
#a_required_path = join(parent_path, 'lib/Clinicasurfstat')
input_directory = join(split(realpath(__file__))[0], 'external-data/clinica_surfstat')
csv_file  = join(split(realpath(__file__))[0], 'external-data/clinica_surfstat/csv_file/template.csv')
str_format = '%s %s %s %f'
linear_model = '1 + Label + Gender + Age'
#output_directory = tempfile.mkdtemp()
output_directory = "~/test"
print 'Output dir %s' % output_directory
contrast = 'Label'
start = time.time()
surfstat = clinica_surfstat(input_directory,output_directory, linear_model, contrast, csv_file, str_format)
surfstat.run()
time_consuming = time.time() - start
print 'END! time consuming is : %s' %time_consuming
