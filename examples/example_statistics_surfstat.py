# -*- coding: utf-8 -*-
"""
=================================================
Example of cilinca_surfstat 
=================================================

This example is to elaborate how to do some statistical analysis(GLM) for the preprocessed data, here, our null hypothesis is :
The thickness between AD and CN is the same.

the parameter explanation is below:
      :param: linear_model: string, the linear model that fit into the GLM, for example '1 + Label + Gender + Age'.
      :param: str_format: string, the format that you want to use for your CSV file column variables, it depends on your CSV file.
      :param: tsv_file: string, the path to your csv file.
      :param: caps_directory:  the output file from recon-all pipeline,specifically, files: ?h.thickness.fwhm**.mgh.
              we put in the directory 'data/Recon-all_Output'.
      :param: contrast:  string, depending on what you want to do, there are two kinds of contrast, one is categorized facor contrast, like 'Label',
              this will return you 8 images, including positive contrast results and negative contrast result; another is continuous factor result, 
              like 'age', which will return you 4 images.
      :param: output_directory: the directory to contain the result images. 
       Defaut parameters, we set these parameters to be some default values, but you can also set it by yourself:
      :param: size_o
      f_fwhm: fwhm for the surface smoothing, default is 20, integer.
      :param: threshold_uncorrected_pvalue: threshold to display the uncorrected Pvalue, float.
      :param: threshold_corrected_pvalue: the threshold to display the corrected cluster, default is 0.05, float.
      :param: cluster_threshold: threshold to define a cluster in the process of cluster-wise correction, default is 0.001, float.
      
      
Outputs:
      after the clinica_surfstat pipeline, we will get the results images in the output_directory, also, in the output_directory, we will also
      have a log file 'matlab_output.log', which includes the matlab version information and the surfstat progress information.

Note: as we will use OpenGL to render the result images, and after Matlab2014, they changed the opengl algorithms to make rendering more flexible, 
      meanwhile, maybe a little slower than the older version(not always), and we always recommend using the hardware for OpenGL, which is default
      mode in clinica_surfstat. If you have more than more matlab version in your system, to choose which matlab version that you want to use in your local machine, you should export an environment variable
      'MATLABCMD' in your bashrc file to point to the needed matlab version, if 'MATLABCMD' is not defined, clinica_surfstat will use default matlab
      command line 'matlab'.
      For Mac os x, opengl software mode is not supported, so it will always be opengl hardware mode.

@author: Junhao WEN
"""

from __future__ import absolute_import
from clinica.pipeline.statistics.surfstat import clinica_surfstat
import tempfile
import time
# from os.path import realpath, split, join

# TODO adapt my code to non-BIDS version , like in clinica-data

caps_dir =
tsv_file  =
str_format =
linear_model =

print 'Output dir is in the same CAPS folder with the input'
contrast = 'group'
group_label = 'test'
analysis_series_id = 'default'
working_directory='~/test'
start = time.time()
surfstat = clinica_surfstat(caps_dir, tsv_file, linear_model, contrast, str_format, group_label, working_directory=working_directory)
surfstat.run("MultiProc", plugin_args={'n_procs': 4})
time_consuming = time.time() - start
print 'END! time consuming is : %s' %time_consuming

# command line example:
# clinica run statistics-surfstat /Volumes/dataARAMIS/users/CLINICA/CLINICA_datasets/for_testing/test_surfstat /Volumes/dataARAMIS/users/CLINICA/CLINICA_datasets/for_testing/test_surfstat/analysis-series-default/subjects/subjects_group_list.tsv '1 + group + sex + age' 'group' '%s %s %s %f' 'test'

