#!/usr/bin/python#
# -*- coding: utf-8 -*-
"""
Created on 12/08/2016

@author: Pietro Gori
"""
from __future__ import absolute_import

def weighted_prototypes(input_directory, output_directory, linear_model, contrast, csv_file, str_format):
    """
        This is to use surfstat to do the Group analysis for the reconAll outputs, after the reconAll pipeline, you should just define the paths to
        surfstatGroupAnalysis, and create the CSV file, and run the pipeline, at last, you will get the results images.

        Inputnode
        ---------
        surfstat
        Inputs: input_directory:  the output file from recon-all pipeline,specifically, files: ?h.thickness.fwhm**.mgh.
                output_directory: the directory to contain the result images.
                linear_model: string, the linear model that fit into the GLM, for example '1+Lable'.
                contrast: string, the contrast matrix for GLM, if the factor you choose is categorized variable, clinica_surfstat will create two contrasts,
                          for example, contrast = 'Label', this will create contrastpos = Label.AD - Label.CN, contrastneg = Label.CN - Label.AD; if the fac-
                          tory that you choose is a continuous factor, clinica_surfstat will just create one contrast, for example, contrast = 'Age', but note,
                          the string name that you choose should be exactly the same with the columns names in your csv_file.
                csv_file: string, the path to your csv file.
                str_format: string, the str_format which uses to read your csv file, the typy of the string should corresponds exactly with the columns in the csv file.
                path_to_matscript: this is the path to find the matlabscript where we put it in clinica/lib/clinicasurfstat

          For more infomation about SurfStat, please check:
          http://www.math.mcgill.ca/keith/surfstat/

        Outputnode:

        :param: input_directory:  the output file from recon-all pipeline,specifically, files: ?h.thickness.fwhm**.mgh.
        :param: output_directory: the directory to contain the result images.
        :param: linear_model: string, the linear model that fit into the GLM, for example '1+Lable'.
        :param: contrast: string, the contrast matrix for GLM, if the factor you choose is categorized variable, clinica_surfstat will create two contrasts,
                          for example, contrast = 'Label', this will create contrastpos = Label.AD - Label.CN, contrastneg = Label.CN - Label.AD; if the fac-
                          tory that you choose is a continuous factor, clinica_surfstat will just create one contrast, for example, contrast = 'Age', but note,
                          the string name that you choose should be exactly the same with the columns names in your csv_file.
        :param: csv_file: string, the path to your csv file.
        :param: str_format: string, the str_format which uses to read your csv file, the typy of the string should corresponds exactly with the columns in the csv file.

        return images in output_directory

    """

    from nipype.interfaces.utility import Function
    import nipype.pipeline.engine as pe
    from os.path import realpath, split, join, dirname

    cwd_path = split(realpath(__file__))[0]
    parent_path = dirname(dirname(cwd_path))
    path_to_matscript = join(parent_path, 'lib/clinicasurfstat')

    def runmatlab(input_directory, output_directory, linear_model, contrast, csv_file, str_format, path_to_matscript):
        from nipype.interfaces.matlab import MatlabCommand, get_matlab_command
        from os.path import join

        MatlabCommand.set_default_matlab_cmd(get_matlab_command())#this is to set the matlab_path(os.environ) in your bashrc file, to choose which version of matlab do you wanna use
        matlab = MatlabCommand()

        # add the dynamic traits
        #openGL_trait = traits.Bool(True, argstr='-nosoftwareopengl', usedefault=True, desc='Switch on hardware openGL', nohash=True)
        #matlab.input_spec.add_trait(matlab.input_spec(), 'nosoftwareopengl', openGL_trait() )

        matlab.inputs.args = '-nosoftwareopengl'
        matlab.inputs.paths = path_to_matscript  #CLINICA_HOME, this is the path to add into matlab, addpath
        matlab.inputs.script = """
        clinicasurfstat('%s', '%s', '%s', '%s', '%s', '%s');
        """%(input_directory, output_directory, linear_model, contrast, csv_file, str_format)  # here, we should define the inputs for the matlab function that you want to use
        matlab.inputs.mfile = True # this will create a file: pyscript.m , the pyscript.m is the default name
        matlab.inputs.single_comp_thread = False  #this will stop runing with single thread
        matlab.inputs.logfile = join(output_directory, "matlab_output.log")
        print "matlab logfile is located in the folder: %s" % matlab.inputs.logfile
        print "matlab script command = %s" % matlab.inputs.script
        print "MatlabCommand inputs flag: single_comp_thread = %s" % matlab.inputs.single_comp_thread
        print "MatlabCommand choose which matlab to use: %s" % get_matlab_command()
        print "MatlabCommand inputs flag: nosoftwareopengl = %s" % matlab.inputs.args
        out = matlab.run()
        return out

    surfstat = pe.Node(name='surfstat',
                   interface=Function(input_names=['input_directory', 'output_directory', 'linear_model',
                                         'contrast', 'csv_file', 'str_format', 'path_to_matscript'],
                                      output_names=[ ],
                                      function=runmatlab))
    surfstat.inputs.input_directory = input_directory
    surfstat.inputs.output_directory = output_directory
    surfstat.inputs.linear_model = linear_model
    surfstat.inputs.contrast = contrast
    surfstat.inputs.csv_file = csv_file
    surfstat.inputs.str_format = str_format
    surfstat.inputs.path_to_matscript = path_to_matscript

    return surfstat
