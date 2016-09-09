#!/usr/bin/python#
# -*- coding: utf-8 -*-
"""
Created on 12/08/2016

=================================================
Example of weighted_prototypes.py
=================================================

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
bound_limit_input=0
degree_precision_input=0
num_iter_modularity_input=0
minimum_number_fibers_cluster_input=1
minValueTau_input=0
increase_radius_input=0

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
