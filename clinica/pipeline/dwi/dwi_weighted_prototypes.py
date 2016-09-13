#!/usr/bin/python#
# -*- coding: utf-8 -*-
"""
Created on 12/08/2016

% Compute Streamline Weighted Prototypes of a bundle in VTK format
%
% Usage: weighted_prototypes(working_dir,filename_bundle,lambda_g,lambda_a,lambda_b,path_matlab_functions,path_CPP_code,path_Community_latest,bound_limit_input,degree_precision_input,num_iter_modularity_input,minimum_number_fibers_cluster_input,minValueTau_input,increase_radius_input)
%
% MANDATORY INPUTS:
% - clinica_path: absolute path to clinica home
% - working_dir: directory where files are saved
% - filename_bundle: filename of the fiber bundle which must be a .vtk file.
% It must have a list of 3D points and then it should use the keyword
% "LINES" for polygons. Each row describes a streamline. The first number
% is the number of points. The other numbers are the indexes of the points
% previously listed.
% - lambda_g: geometric kernel bandidth (as in usual currents)
% - lambda_a: kernel bandwidth of the STARTING structure
% - lambda_b: kernel bandwidth of the ENDING structure
% - path_matlab_functions
% - path_cpp_code: absolute path to C++ code
% - path_community_latest: absolute path to Community detection folder
% To note: streamlines must have a consistent orientation, namely they must
% have the same starting and ending ROIs
%
% OPTIONAL INPUTS:
% - bound_limit_input: maximum average angle (in radians) that a streamline
% may have with the other streamlines in the framework of weighted currents.
% Default value is 1.5359 = 88 degrees
% - degree_precision_input: percentage of the norm of the bundle explained by the
% weighted prototypes. Default value is 0.15, which means that the weighted
% prototypes will explain (1-0.15)*100 % of the norm of the bundle in the
% framework of weighted currents.
% - num_iter_modularity_input: Modularity computation is based on a greedy
% approach. Results may differ between iterations. The greater number of
% iterations, the better. Default value is 10
% See "Fast unfolding of community hier archies in large networks", V. Blondel et al.
% - minimum_number_fibers_cluster_input: Clustering based on modularity may
% result in unbalanced clusters. We remove the clusters which have less
% than minimum_number_fibers_cluster_input fibers. Default value is 10
% - minValueTau_input: We remove the prototypes that approximate less than
% minValueTau_input fibers. Default value is 1
% - increase_radius_input: All tubes are normalised such that the maximum
% radius is equal to 1mm. We then augment all radii of
% increase_radius_input. Default value is 0.02
%
% This function requires:
% - The binary files of the C++ functions in the folder cpp_code
% - CMake > 2.8
% - VTK > 6
% - ITK
% - Louvain community detection (https://sites.google.com/site/findcommunities/newversion/community.tgz?attredirects=0)
% which is already present
% - Eigen (http://eigen.tuxfamily.org/index.php?title=Main_Page)
% which is also already present

@author: pietro.gori
"""
from __future__ import absolute_import

def weighted_prototypes(clinica_path,working_dir,filename_bundle,lambda_g,lambda_a,lambda_b,bound_limit_input,degree_precision_input,num_iter_modularity_input,minimum_number_fibers_cluster_input,minValueTau_input,increase_radius_input):

    from nipype.interfaces.utility import Function
    import nipype.pipeline.engine as pe
    from os.path import realpath, split, join, dirname

    path_to_matscript = join(clinica_path, 'clinica/lib/weighted_prototypes_lib') # concatenate parent directory path and 'lib/WeightedPrototypes'
    path_matlab_functions = join(path_to_matscript, 'matlab_functions')
    path_cpp_code = join(path_to_matscript, 'cpp_code')
    path_community_latest = join(path_to_matscript, 'community_latest')

    # DEBUG
    print 'Folder of weighted_prototypes_lib is : %s' % path_to_matscript
    print 'Folder of matlab_functions is : %s' % path_matlab_functions
    print 'Folder of cpp_code is : %s' % path_cpp_code
    print 'Folder of community_latest is : %s' % path_community_latest
      
    def runmatlab(path_to_matscript,working_dir,filename_bundle,lambda_g,lambda_a,lambda_b,path_matlab_functions,path_cpp_code,path_community_latest,bound_limit_input,degree_precision_input,num_iter_modularity_input,minimum_number_fibers_cluster_input,minValueTau_input,increase_radius_input):
                      
        from nipype.interfaces.matlab import MatlabCommand, get_matlab_command
        from os.path import join
        import sys, os
        
        # here, we check out the os, basically, clinica works for linux and MAC OS X.
        if sys.platform.startswith('linux'):
            print "###Note: your platform is linux, the default command line for Matlab(matlab_cmd) is matlab, but you can also export a variable MATLABCMD,  which points to your matlab,  in your .bashrc to set matlab_cmd, this can help you to choose which Matlab to run when you have more than one Matlab. "
        elif sys.platform.startswith('darwin'):
            try:
                if 'MATLABCMD' not in os.environ:
                    raise RuntimeError("###Note: your platform is MAC OS X, the default command line for Matlab(matlab_cmd) is matlab, but it does not work on OS X, you mush export a variable MATLABCMD, which points to your matlab, in your .bashrc to set matlab_cmd.")
            except Exception as e:
                print(str(e))
                exit(1)            
        else:
            print "Clinica will not work on your platform "
            
        MatlabCommand.set_default_matlab_cmd(get_matlab_command())#this is to set the matlab_path(os.environ) in your bashrc file, to choose which version of matlab do you wanna use     
        matlab = MatlabCommand()
        matlab.inputs.args = '-nosoftwareopengl'
        matlab.inputs.paths = path_to_matscript
        matlab.inputs.script = """weighted_prototypes('%s','%s', %f, %f, %f, '%s', '%s', '%s', %f, %f, %d, %d, %f, %f);"""%(working_dir,filename_bundle,lambda_g,lambda_a,lambda_b,path_matlab_functions,path_cpp_code,path_community_latest,bound_limit_input,degree_precision_input,num_iter_modularity_input,minimum_number_fibers_cluster_input,minValueTau_input,increase_radius_input)
        matlab.inputs.mfile = True # this will create a file: pyscript.m , the pyscript.m is the default name
        matlab.inputs.single_comp_thread = False  #this will stop runing with single thread  
        matlab.inputs.logfile = join(working_dir, "matlab_output.log")
        print "matlab logfile is located in : %s" % matlab.inputs.logfile
        print "matlab script command = %s" % matlab.inputs.script
        print "MatlabCommand inputs flag: single_comp_thread = %s" % matlab.inputs.single_comp_thread
        print "MatlabCommand choose which matlab to use(matlab_cmd): %s" % get_matlab_command()
        print "MatlabCommand inputs flag: nosoftwareopengl = %s" % matlab.inputs.args

        out = matlab.run()
        return out
    # end def runmatlab
        
    WP = pe.Node(name='WP', interface=Function(input_names=['path_to_matscript','working_dir','filename_bundle','lambda_g','lambda_a','lambda_b','path_matlab_functions','path_cpp_code','path_community_latest','bound_limit_input','degree_precision_input','num_iter_modularity_input','minimum_number_fibers_cluster_input','minValueTau_input','increase_radius_input'], output_names=[ ], function=runmatlab))
    WP.inputs.path_to_matscript = path_to_matscript
    WP.inputs.working_dir = working_dir
    WP.inputs.filename_bundle = filename_bundle
    WP.inputs.lambda_g = lambda_g
    WP.inputs.lambda_a = lambda_a    
    WP.inputs.lambda_b = lambda_b
    WP.inputs.path_matlab_functions = path_matlab_functions
    WP.inputs.path_cpp_code = path_cpp_code
    WP.inputs.path_community_latest = path_community_latest
    WP.inputs.bound_limit_input = bound_limit_input
    WP.inputs.degree_precision_input = degree_precision_input
    WP.inputs.num_iter_modularity_input = num_iter_modularity_input
    WP.inputs.minimum_number_fibers_cluster_input = minimum_number_fibers_cluster_input
    WP.inputs.minValueTau_input = minValueTau_input
    WP.inputs.increase_radius_input = increase_radius_input
    
    return WP