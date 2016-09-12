#!/usr/bin/python#
# -*- coding: utf-8 -*-
"""
Created on 12/08/2016

=================================================
Example of dwi_weighted_prototypes.py
=================================================

This example is to elaborate how to do some statistical analysis(GLM) for the preprocessed data, here, our null hypothesis is :
The thickness between AD and CN is the same.

the parameter explanation is below:
      :param: linear_model: string, the linear model that fit into the GLM, for example '1 + Label + Gender + Age'.
      :param: str_format: string, the format that you want to use for your CSV file column variables, it depends on your CSV file.
      :param: csv_file: string, the path to your csv file.
      :param: input_directory:  the output file from recon-all pipeline,specifically, files: ?h.thickness.fwhm**.mgh.
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


@author: pietro.gori

"""

from __future__ import absolute_import
from clinica.pipeline.dwi.dwi_weighted_prototypes import weighted_prototypes
from os.path import realpath, split, join
from os import makedirs
import tempfile
import time
import errno

cwd_path = split(realpath(__file__))[0]

## MANDATORY PARAMETERS ##
lambda_g=4.0
lambda_a=3.0
lambda_b=3.0
type='small' # It can be 'small' 'medium' or 'big'
filename_bundle=join(cwd_path, 'external-data/dwi_weighted_prototypes/bundle_%s.vtk') % type

## OPTIONAL PARAMETERS ##
bound_limit_input=0.0
degree_precision_input=0.0
num_iter_modularity_input=0
minimum_number_fibers_cluster_input=1
minValueTau_input=0.0
increase_radius_input=0.0

working_dir = join(cwd_path, 'external-data/dwi_weighted_prototypes/example_bundle_%s_lambda_g_%.1f_lambda_a_%.1f_lambda_b_%.1f') % (type,lambda_g,lambda_a,lambda_b)
try:
    makedirs(working_dir)
except OSError as exc:
    if exc.errno == errno.EEXIST:
        print 'WARNING: %s is already present. It seems that you have already run this experiment. Delete the folder to proceed.' % working_dir
        raise exc
    pass

print 'Working directory: %s' % working_dir

start = time.time()
WP = weighted_prototypes(working_dir,filename_bundle,lambda_g,lambda_a,lambda_b,bound_limit_input,degree_precision_input,num_iter_modularity_input,minimum_number_fibers_cluster_input,minValueTau_input,increase_radius_input)
WP.run()
time_consuming = time.time() - start

print 'END! Weighted Prototypes computed in : %s' %time_consuming
