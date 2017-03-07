#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""This module contains functions used for the surfstat() pipeline"""
import os

__author__ = "Junhao Wen, Alexandre Routier"
__copyright__ = "Copyright 2016, The Aramis Lab Team"
__credits__ = ["Michael Bacci", "Junhao Wen"]
__license__ = "??"
__version__ = "1.0.0"
__maintainer__ = "Junhao Wen"
__email__ = "junhao.Wen@inria.fr"
__status__ = "Development"

# transfer any path to be absolute path.
def absolute_path(arg):
    """
    Fetch the absolute path for any input path

    :param arg:
    :return:
    """
    if arg[:1] == '~':
        return os.path.expanduser(arg)
    elif arg[:1] == '.':
        return os.getcwd()
    else:
        return os.path.join(os.getcwd(), arg)


def get_vars(input_directory, subjects_visits_tsv, group_label):
    """
    Fetch all the intermedial variables for this workflow

    :param input_directory:
    :param subjects_visits_tsv:
    :param group_label:
    :return:
    """
    import os
    from glob import glob
    from shutil import copy
    import clinica.pipeline as clp

    path_to_matscript = os.path.join(os.path.dirname(clp.__path__[0]), 'lib/clinicasurfstat')

    # CAPS input and output vars
    input_directory = os.path.expanduser(input_directory)
    surfstat_input_dir = os.path.join(input_directory, 'subjects')

    group_id = 'group-' + group_label
    output_directory = os.path.join(input_directory, 'groups', group_id, 'statistics/surfstat/clinica-surfstat')
    if not os.path.exists(output_directory):
        try:
            os.makedirs(output_directory)
        except:
            raise OSError("Surfstat: can't create destination directory (%s)!" % (output_directory))

    # cp the subjects_visits_tsv to the result folder
    copied_tsv = output_directory + '/subjects_group_list.tsv'
    copy(subjects_visits_tsv, copied_tsv)

    return path_to_matscript, surfstat_input_dir, output_directory

def runmatlab(input_directory,
              output_directory,
              subjects_visits_tsv,
              design_matrix, contrast,
              str_format,
              glm_type,
              path_to_matscript,
              full_width_at_half_maximum,
              threshold_uncorrected_pvalue,
              threshold_corrected_pvalue,
              cluster_threshold):
    """
    a wrapper the matlab script of surfstat with nipype.

    :param input_directory:
    :param output_directory:
    :param subjects_visits_tsv:
    :param design_matrix:
    :param contrast:
    :param str_format:
    :param path_to_matscript:
    :param full_width_at_half_maximum:
    :param threshold_uncorrected_pvalue:
    :param threshold_corrected_pvalue:
    :param cluster_threshold:
    :return:
    """
    from nipype.interfaces.matlab import MatlabCommand, get_matlab_command
    from os.path import join
    import sys, os
    # here, we check out the os, basically, clinica works for linux and MAC OS X.
    if sys.platform.startswith('linux'):
        print "###Note: your platform is linux, the default command line for Matlab(matlab_cmd) is matlab, but you can also export a variable MATLABCMD,  which points to your matlab,  in your .bashrc to set matlab_cmd, this can help you to choose which Matlab to run when you have more than one Matlab. "
    elif sys.platform.startswith('darwin'):
        try:
            if not 'MATLABCMD' in os.environ:
                raise RuntimeError(
                    "###Note: your platform is MAC OS X, the default command line for Matlab(matlab_cmd) is matlab, but it does not work on OS X, you mush export a variable MATLABCMD, which points to your matlab, in your .bashrc to set matlab_cmd. Note, Mac os x will always choose to use OpengGl hardware mode.")
        except Exception as e:
            print(str(e))
            exit(1)
    else:
        print "Clinica will not work on your platform "

    MatlabCommand.set_default_matlab_cmd(
        get_matlab_command())  # this is to set the matlab_path(os.environ) in your bashrc file, to choose which version of matlab do you wanna use
    # here, set_default_matlab_cmd is a @classmethod
    matlab = MatlabCommand()

    # add the dynamic traits
    # openGL_trait = traits.Bool(True, argstr='-nosoftwareopengl', usedefault=True, desc='Switch on hardware openGL', nohash=True)
    # matlab.input_spec.add_trait(matlab.input_spec(), 'nosoftwareopengl', openGL_trait() )
    if sys.platform.startswith('linux'):
        matlab.inputs.args = '-nosoftwareopengl'  # Bug, for my laptop, it does not work, but the command line does have the flag -nosoftwareopengl, we should try on other computer's matlab to check if this flag works!
    matlab.inputs.paths = path_to_matscript  # CLINICA_HOME, this is the path to add into matlab, addpath

    matlab.inputs.script = """
    clinicasurfstat('%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', %d, '%s', %.3f, '%s', %.3f, '%s', %.3f);
    """ % (input_directory, output_directory, subjects_visits_tsv, design_matrix, contrast, str_format, glm_type, 'sizeoffwhm',
           full_width_at_half_maximum,
           'thresholduncorrectedpvalue', threshold_uncorrected_pvalue, 'thresholdcorrectedpvalue',
           threshold_corrected_pvalue, 'clusterthreshold',
           cluster_threshold)  # here, we should define the inputs for the matlab function that you want to use
    matlab.inputs.mfile = True  # this will create a file: pyscript.m , the pyscript.m is the default name
    matlab.inputs.single_comp_thread = False  # this will stop runing with single thread
    matlab.inputs.logfile = join(output_directory, "matlab_output.log")
    print "Matlab logfile is located in the folder: %s" % matlab.inputs.logfile
    print "Matlab script command = %s" % matlab.inputs.script
    print "MatlabCommand inputs flag: single_comp_thread = %s" % matlab.inputs.single_comp_thread
    print "MatlabCommand choose which matlab to use(matlab_cmd): %s" % get_matlab_command()
    if sys.platform.startswith('linux'):
        print "MatlabCommand inputs flag: nosoftwareopengl = %s" % matlab.inputs.args

    out = matlab.run()
    return out