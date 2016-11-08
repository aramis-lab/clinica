#!/usr/bin/python#
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 15:20:40 2016

@author: Junhao WEN
"""
from __future__ import absolute_import

def clinica_surfstat(input_directory,
                     subjects_visits_tsv,
                     linear_model,
                     contrast,
                     str_format,
                     group_label,
                     size_of_fwhm = 20,
                     threshold_uncorrected_pvalue = 0.001,
                     threshold_corrected_pvalue = 0.050,
                     cluster_threshold = 0.001,
                     working_directory=None):
    """
        This is to use surfstat to do the Group analysis for the reconAll outputs, after the reconAll pipeline, you should just define the paths to
        surfstatGroupAnalysis, and create the tsv file, and run the pipeline, at last, you will get the results images.

        Inputs
        ---------
        surfstat
        Inputs: :param input_directory:  the output file from recon-all pipeline,specifically, files: ?h.thickness.fwhm**.mgh.
                :param linear_model: string, the linear model that fits into the GLM, for example '1+Lable'.
                :param contrast: string, the contrast matrix for GLM, if the factor you choose is categorized variable, clinica_surfstat will create two contrasts,
                          for example, contrast = 'Label', this will create contrastpos = Label.AD - Label.CN, contrastneg = Label.CN - Label.AD; if the fac-
                          tory that you choose is a continuous factor, clinica_surfstat will just create one contrast, for example, contrast = 'Age', but note,
                          the string name that you choose should be exactly the same with the columns names in your subjects_visits_tsv.
                :param subjects_visits_tsv: string, the path to your tsv file.
                :param str_format: string, the str_format which uses to read your tsv file, the typy of the string should corresponds exactly with the columns in the tsv file.
                 Defaut parameters, we set these parameters to be some default values, but you can also set it by yourself:
                :param size_of_fwhm: fwhm for the surface smoothing, default is 20, integer.
                :param threshold_uncorrected_pvalue: threshold to display the uncorrected Pvalue, float, default is 0.001.
                :param threshold_corrected_pvalue: the threshold to display the corrected cluster, default is 0.05, float.
                :param cluster_threshold: threshold to define a cluster in the process of cluster-wise correction, default is 0.001, float.
                :param: working_directory: define where to put the infomation of the nipype workflow.

          For more infomation about SurfStat, please check:
          http://www.math.mcgill.ca/keith/surfstat/

        Outputs:
        return result images in output_directory of clinicasurfstat matlab script.

    """
    
    from nipype.interfaces.utility import Function
    import nipype.pipeline.engine as pe
    from glob import glob
    import os
    import nipype.interfaces.utility as niu
    from shutil import copy
    from tempfile import mkdtemp
    
    cwd_path = os.path.split(os.path.realpath(__file__))[0]
    parent_path = os.path.dirname(os.path.dirname(cwd_path))
    path_to_matscript = os.path.join(parent_path, 'lib/clinicasurfstat')

    #transfer any path to be absolute path.
    def absolute_path(arg):
        if arg[:1] == '~':
            return os.path.expanduser(arg)
        elif arg[:1] == '.':
            return os.getcwd()
        else:
            return os.path.join(os.getcwd(), arg)

    def CAPS_input(input_directory):
        # from glob import glob
        # import os
        input_inter_path = glob(os.path.join(input_directory, '*/subjects'))
        return input_inter_path[0]

    def CAPS_output(input_directory):
        intermediate_path = glob(os.path.join(input_directory, '*'))
        # analysis_id = 'group-' + str(size_of_fwhm) + '-' + str(threshold_uncorrected_pvalue) + '-' + str(threshold_corrected_pvalue) + '-' + str(cluster_threshold)
        analysis_id = 'group-' + group_label
        output_inter_path = os.path.join(intermediate_path[0], 'group', analysis_id, 'statistics/surfstat/clinica-surfstat')
        if not os.path.exists(output_inter_path):
            try:
                os.makedirs(output_inter_path)
            except:
                raise OSError("Surfstat: can't create destination directory (%s)!" % (output_inter_path))
        return output_inter_path

    output_directory= CAPS_output(input_directory)

    def runmatlab(input_directory, output_directory, subjects_visits_tsv, linear_model, contrast, str_format, path_to_matscript,
                  size_of_fwhm, threshold_uncorrected_pvalue, threshold_corrected_pvalue, cluster_threshold ):
        from nipype.interfaces.matlab import MatlabCommand, get_matlab_command
        from os.path import join
        import sys, os
        # here, we check out the os, basically, clinica works for linux and MAC OS X.
        if sys.platform.startswith('linux'):
            print "###Note: your platform is linux, the default command line for Matlab(matlab_cmd) is matlab, but you can also export a variable MATLABCMD,  which points to your matlab,  in your .bashrc to set matlab_cmd, this can help you to choose which Matlab to run when you have more than one Matlab. "
        elif sys.platform.startswith('darwin'):
            try:
                if not 'MATLABCMD' in  os.environ:
                    raise RuntimeError("###Note: your platform is MAC OS X, the default command line for Matlab(matlab_cmd) is matlab, but it does not work on OS X, you mush export a variable MATLABCMD, which points to your matlab, in your .bashrc to set matlab_cmd. Note, Mac os x will always choose to use OpengGl hardware mode.")
            except Exception as e:
                print(str(e))
                exit(1)            
        else:
            print "Clinica will not work on your platform "
            
        MatlabCommand.set_default_matlab_cmd(get_matlab_command())#this is to set the matlab_path(os.environ) in your bashrc file, to choose which version of matlab do you wanna use     
        # here, set_default_matlab_cmd is a @classmethod
        matlab = MatlabCommand()
        
        # add the dynamic traits
        #openGL_trait = traits.Bool(True, argstr='-nosoftwareopengl', usedefault=True, desc='Switch on hardware openGL', nohash=True)
        #matlab.input_spec.add_trait(matlab.input_spec(), 'nosoftwareopengl', openGL_trait() )
        if sys.platform.startswith('linux'):
            matlab.inputs.args = '-nosoftwareopengl' # Bug, for my laptop, it does not work, but the command line does have the flag -nosoftwareopengl, we should try on other computer's matlab to check if this flag works!
        matlab.inputs.paths = path_to_matscript  #CLINICA_HOME, this is the path to add into matlab, addpath
        # there is a bug to transfer the diff type of inputs througn the interface to matlab script, they always make it to be a string, that is because .script is a Str trait. That is because here, we use the nipype existed intreface, if you write your own interface, I think you can choose the type of your 
        # variables that you want to transfer to the matlab script.
        matlab.inputs.script = """
        clinicasurfstat('%s', '%s', '%s', '%s', '%s', '%s', '%s', %d, '%s', %.3f, '%s', %.3f, '%s', %.3f);
        """%(input_directory, output_directory, subjects_visits_tsv, linear_model, contrast, str_format, 'sizeoffwhm', size_of_fwhm,
             'thresholduncorrectedpvalue', threshold_uncorrected_pvalue, 'thresholdcorrectedpvalue', threshold_corrected_pvalue, 'clusterthreshold', cluster_threshold)  # here, we should define the inputs for the matlab function that you want to use
        matlab.inputs.mfile = True # this will create a file: pyscript.m , the pyscript.m is the default name
        matlab.inputs.single_comp_thread = False  #this will stop runing with single thread  
        matlab.inputs.logfile = join(output_directory, "matlab_output.log")
        print "Matlab logfile is located in the folder: %s" % matlab.inputs.logfile
        print "Matlab script command = %s" % matlab.inputs.script
        print "MatlabCommand inputs flag: single_comp_thread = %s" % matlab.inputs.single_comp_thread
        print "MatlabCommand choose which matlab to use(matlab_cmd): %s" % get_matlab_command()
        if sys.platform.startswith('linux'):        
            print "MatlabCommand inputs flag: nosoftwareopengl = %s" % matlab.inputs.args

        out = matlab.run()
        return out

    outputnode = pe.Node(niu.IdentityInterface(fields=['surfstat_result']), name='outputnode')
    surfstat = pe.Node(name='surfstat',
                   interface=Function(input_names=['input_directory', 'output_directory', 'subjects_visits_tsv', 'linear_model',
                                         'contrast', 'str_format', 'path_to_matscript', 'size_of_fwhm', 'threshold_uncorrected_pvalue',
                                         'threshold_corrected_pvalue', 'cluster_threshold'],
                                      output_names=['out_images'],
                                      function=runmatlab))
    surfstat.inputs.input_directory = CAPS_input(input_directory)
    surfstat.inputs.output_directory = output_directory
    surfstat.inputs.linear_model = linear_model    
    surfstat.inputs.contrast = contrast
    surfstat.inputs.subjects_visits_tsv = subjects_visits_tsv
    surfstat.inputs.str_format = str_format
    surfstat.inputs.path_to_matscript = path_to_matscript
    surfstat.inputs.size_of_fwhm = size_of_fwhm
    surfstat.inputs.threshold_uncorrected_pvalue = threshold_uncorrected_pvalue
    surfstat.inputs.threshold_corrected_pvalue = threshold_corrected_pvalue
    surfstat.inputs.cluster_threshold = cluster_threshold

    # cp the subjects_visits_tsv to the result folder
    copied_tsv = output_directory + '/subjects_group_list.tsv'
    copy(subjects_visits_tsv, copied_tsv)

    if working_directory is None:
        working_directory = mkdtemp()
    else:
        working_directory = absolute_path(working_directory)

    surfstat_wf = pe.Workflow(name='surfstat_workflow', base_dir=working_directory)
    surfstat_wf.connect(surfstat, 'out_images', outputnode, 'surfstat_result')

    return surfstat_wf
