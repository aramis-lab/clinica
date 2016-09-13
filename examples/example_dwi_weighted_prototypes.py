#!/usr/bin/python#
# -*- coding: utf-8 -*-
"""
Created on 12/08/2016

=================================================
Example of dwi_weighted_prototypes.py
=================================================

Example about the approximation of a streamline bundle into weighted prototypes.

 MANDATORY INPUTS:
 - lambda_g: geometric kernel bandidth (as in usual currents)
 - lambda_a: kernel bandwidth of the STARTING structure
 - lambda_b: kernel bandwidth of the ENDING structure
 - type: Type of bundle to approximate. It can be 'small' 'medium' or 'big'

 OPTIONAL INPUTS:
 - bound_limit_input: maximum average angle (in radians) that a streamline may have with the other streamlines in the
                       framework of weighted currents. Default value is 1.5359 = 88 degrees
 - degree_precision_input: percentage of the norm of the bundle explained by the weighted prototypes.
                            Default value is 0.15, which means that the weighted prototypes will explain (1-0.15)*100
                            of the norm of the bundle in the framework of weighted currents.
 - num_iter_modularity_input: Modularity computation is based on a greedy approach. Results may differ between
                               iterations. The greater number of iterations, the better. Default value is 10
                               See "Fast unfolding of community hier archies in large networks", V. Blondel et al.
 - minimum_number_fibers_cluster_input: Clustering based on modularity may result in unbalanced clusters.
                                         We remove the clusters which have less than minimum_number_fibers_cluster_input
                                         fibers. Default value is 10
 - minValueTau_input: We remove the prototypes that approximate less than minValueTau_input fibers. Default value is 1
 - increase_radius_input: All tubes are normalised such that the maximum radius is equal to 1mm. We then augment all
                           radii of increase_radius_input. Default value is 0.02

Outputs:
    - graph.* : Binary files containing the Gramiam
    - Clusters.vtk: VTK file containing the original bundle divided in fascicles (clusters)
    - NoOutlier.vtk: VTK file containing the original bundle without the streamlines considered as outliers
    - Prototypes.vtk: VTK files containing the weighted prototypes
    - Prototypes_tubes.vtk: VTK files containing the weighted prototypes represented as tubes
    - *.log: Log files containing the outputs of the different steps

This function requires:
 - The binary files of the C++ functions in the folder cpp_code/bin
 - CMake > 2.8
 - VTK > 6
 - ITK
 - Louvain community detection (https://sites.google.com/site/findcommunities/newversion/community.tgz?attredirects=0)
   which is already present, the user just needs to do 'make'
 - Eigen (http://eigen.tuxfamily.org/index.php?title=Main_Page)
   which is also already present

@author: pietro.gori

"""

from __future__ import absolute_import
from clinica.pipeline.dwi.dwi_weighted_prototypes import weighted_prototypes
from os.path import realpath, split, join, dirname
from os import makedirs
import time
import errno

## MANDATORY PARAMETERS ##
lambda_g=4.0
lambda_a=3.0
lambda_b=3.0
type='small' # It can be 'small' 'medium' or 'big'
<<<<<<< 4ddf9ffd54b30ff50262476ab7439ac04a12ba09
filename_bundle=join(cwd_path, 'external-data/dwi_weighted_prototypes/bundle_%s.vtk') % type
=======
>>>>>>> Added clinica home path

## OPTIONAL PARAMETERS ##
bound_limit_input=0.0
degree_precision_input=0.0
num_iter_modularity_input=0
minimum_number_fibers_cluster_input=1
minValueTau_input=0.0
increase_radius_input=0.0

# CODE
filename_bundle=join(dirname(__file__), 'external-data/dwi_weighted_prototypes/bundle_%s.vtk') % type
print 'Chosen bundle is : %s' % filename_bundle
clinica_path = dirname(dirname(__file__))
print 'Home of clinica is : %s' % clinica_path

# Working_dir is where the code is run
working_dir = join(dirname(__file__), 'external-data/dwi_weighted_prototypes/example_bundle_%s_lambda_g_%.1f_lambda_a_%.1f_lambda_b_%.1f') % (type,lambda_g,lambda_a,lambda_b)
try:
    makedirs(working_dir)
except OSError as exc:
    if exc.errno == errno.EEXIST:
        print 'WARNING: %s is already present. It seems that you have already run this experiment. Delete the folder to proceed.' % working_dir
        raise exc
    pass

print 'Working directory: %s' % working_dir

start = time.time()
WP = weighted_prototypes(clinica_path,working_dir,filename_bundle,lambda_g,lambda_a,lambda_b,bound_limit_input,degree_precision_input,num_iter_modularity_input,minimum_number_fibers_cluster_input,minValueTau_input,increase_radius_input)
WP.run()
time_consuming = time.time() - start

print 'END! Weighted Prototypes computed in : %s' %time_consuming
