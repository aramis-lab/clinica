#!/usr/bin/python#
# -*- coding: utf-8 -*-
"""
Created on 12/08/2016

This is to use surfstat to do the Group analysis for the reconAll outputs, after the reconAll pipeline, you should just define the paths to
        surfstatGroupAnalysis, and create the CSV file, and run the pipeline, at last, you will get the results images.

        Inputs
        ---------
        surfstat
        Inputs: :param input_directory:  the output file from recon-all pipeline,specifically, files: ?h.thickness.fwhm**.mgh.
                :param output_directory: the directory to contain the result images.
                :param linear_model: string, the linear model that fit into the GLM, for example '1+Lable'.
                :param contrast: string, the contrast matrix for GLM, if the factor you choose is categorized variable, clinica_surfstat will create two contrasts,
                          for example, contrast = 'Label', this will create contrastpos = Label.AD - Label.CN, contrastneg = Label.CN - Label.AD; if the fac-
                          tory that you choose is a continuous factor, clinica_surfstat will just create one contrast, for example, contrast = 'Age', but note,
                          the string name that you choose should be exactly the same with the columns names in your csv_file.
                :param csv_file: string, the path to your csv file.
                :param str_format: string, the str_format which uses to read your csv file, the typy of the string should corresponds exactly with the columns in the csv file.
                 Defaut parameters, we set these parameters to be some default values, but you can also set it by yourself:
                :param size_of_fwhm: fwhm for the surface smoothing, default is 20, integer.
                :param threshold_uncorrected_pvalue: threshold to display the uncorrected Pvalue, float.
                :param threshold_corrected_pvalue: the threshold to display the corrected cluster, default is 0.05, float.
                :param cluster_threshold: threshold to define a cluster in the process of cluster-wise correction, default is 0.001, float.
          For more infomation about SurfStat, please check:
          http://www.math.mcgill.ca/keith/surfstat/

        Outputs:
        return result images in output_directory

@author: pietro.gori
"""
from __future__ import absolute_import

def weighted_prototypes(working_dir,filename_bundle,lambda_g,lambda_a,lambda_b,bound_limit_input,degree_precision_input,num_iter_modularity_input,minimum_number_fibers_cluster_input,minValueTau_input,increase_radius_input):
    """
        TODO

    """

    from nipype.interfaces.utility import Function
    import nipype.pipeline.engine as pe
    from os.path import realpath, split, join, dirname 
    
    cwd_path = split(realpath(__file__))[0] # current working directory path
    parent_path = dirname(dirname(cwd_path)) # cd ../..
    path_to_matscript = join(parent_path, 'lib/weighted_prototypes_lib') # concatenate parent directory path and 'lib/WeightedPrototypes'
    path_matlab_functions = join(path_to_matscript, 'matlab_functions')
    path_cpp_code = join(path_to_matscript, 'cpp_code')
    path_community_latest = join(path_to_matscript, 'community_latest')
      
      
    def runmatlab(path_to_matscript,working_dir,filename_bundle,lambda_g,lambda_a,lambda_b,path_matlab_functions,path_cpp_code,path_community_latest,bound_limit_input,degree_precision_input,num_iter_modularity_input,minimum_number_fibers_cluster_input,minValueTau_input,increase_radius_input):
                      
        from nipype.interfaces.matlab import MatlabCommand, get_matlab_command
        from os.path import join
        import sys, os
        
        # here, we check out the os, basically, clinica works for linux and MAC OS X.
        if sys.platform.startswith('linux'):
            print "###Note: your platform is linux, the default command line for Matlab(matlab_cmd) is matlab, but you can also export a variable MATLABCMD,  which points to your matlab,  in your .bashrc to set matlab_cmd, this can help you to choose which Matlab to run when you have more than one Matlab. "
        elif sys.platform.startswith('darwin'):
            try:
                if not 'MATLABCMD' in  os.environ:
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