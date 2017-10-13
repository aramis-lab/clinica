#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""This module contains functions used for the statistics_surfstat() pipelines"""

__author__ = "Junhao Wen"
__copyright__ = "Copyright 2016, The Aramis Lab Team"
__credits__ = ["Michael Bacci", "Junhao Wen"]
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__maintainer__ = "Junhao Wen"
__email__ = "junhao.Wen@inria.fr"
__status__ = "Development"

def data_prep(input_directory, subjects_visits_tsv, group_label, glm_type):
    """
        Fetch all the intermedial variables for this workflow

    Args:
        input_directory: CAPS directory
        subjects_visits_tsv: tsv file defining the glm model
        group_label: Current group name for this analysis
        glm_type: based on the hypothesis, you should define one of the glm types, "group_comparison", "correlation"

    Returns:

    """
    import os
    from shutil import copy
    import clinica.pipelines as clp
    import sys
    import pandas as pd

    path_to_matscript = os.path.join(os.path.dirname(clp.__path__[0]), 'lib/clinicasurfstat')

    # CAPS input and output vars
    input_directory = os.path.expanduser(input_directory)
    surfstat_input_dir = os.path.join(input_directory, 'subjects')

    group_id = 'group-' + group_label
    statistics_dir_tsv = os.path.join(input_directory, 'groups', group_id, 'statistics', 'participant.tsv')

    if glm_type == "group_comparison":
        output_directory = os.path.join(input_directory, 'groups', group_id, 'statistics', 'surfstat_group_comparison')
        if not os.path.exists(output_directory):
            try:
                os.makedirs(output_directory)
            except:
                raise OSError("Surfstat: can't create destination directory (%s)!" % (output_directory))
    elif glm_type == "correlation":
        output_directory = os.path.join(input_directory, 'groups', group_id, 'statistics', 'surfstat_correlation_analysis')
        if not os.path.exists(output_directory):
            try:
                os.makedirs(output_directory)
            except:
                raise OSError("Surfstat: can't create destination directory (%s)!" % (output_directory))
    else:
        print "The other GLM situations have not been implemented in this pipelines"
        sys.exit()

    # cp the subjects_visits_tsv to the result folder
    # firstly, check if the subjects_visits_tsv has the same info with the participant.tsv in the folder of statistics.
    # if the participant tsv does not exit, cp subjects_visits_tsv in the folder of statistics too, if it is here, compare them.
    if not os.path.isfile(statistics_dir_tsv):
        copy(subjects_visits_tsv, statistics_dir_tsv)
    else:
        ## compare the two tsv files
        participant_df = pd.io.parsers.read_csv(statistics_dir_tsv, sep='\t')
        subjects_visits_tsv_df = pd.io.parsers.read_csv(subjects_visits_tsv, sep='\t')
        participant_list = list(participant_df.participant_id)
        subjects_visits_tsv_list = list(subjects_visits_tsv_df.participant_id)
        dif_list = list(set(participant_list) - set(subjects_visits_tsv_list))
        try:
            len(dif_list) == 0
        except:
            raise ValueError("It seems that this round of analysis does not contain the same subjects where you want to put the results, please check it!")

        group_tsv = 'group-' + group_label + '_participants.tsv'
        copied_tsv = os.path.join(output_directory, group_tsv)
        copy(subjects_visits_tsv, copied_tsv)
    ## point to the path to the json file
    out_json = os.path.join(output_directory, 'group-' + group_label + '_glm.json')

    #### get the FreeSurfer environment variable: FREESURFER_HOME
    freesurfer_home = os.environ["FREESURFER_HOME"]

    return path_to_matscript, surfstat_input_dir, output_directory, freesurfer_home, out_json

def runmatlab(input_directory,
              output_directory,
              subjects_visits_tsv,
              design_matrix, contrast,
              str_format,
              glm_type,
              group_label,
              freesurfer_home,
              surface_file,
              path_to_matscript,
              full_width_at_half_maximum,
              threshold_uncorrected_pvalue,
              threshold_corrected_pvalue,
              cluster_threshold):
    """
        a wrapper the matlab script of surfstat with nipype.

    Args:
        input_directory: surfstat_input_dir where containing all the subjects' output in CAPS directory
        output_directory: output folder to contain the result in CAPS folder
        subjects_visits_tsv: tsv file containing the glm information
        design_matrix: str, the linear model that fits into the GLM, for example '1+group'.
        contrast: string, the contrast matrix for GLM, if the factor you choose is categorized variable, clinica_surfstat will create two contrasts,
                  for example, contrast = 'Label', this will create contrastpos = Label.AD - Label.CN, contrastneg = Label.CN - Label.AD; if the fac-
                  tory that you choose is a continuous factor, clinica_surfstat will just create one contrast, for example, contrast = 'Age', but note,
                  the string name that you choose should be exactly the same with the columns names in your subjects_visits_tsv.
        str_format:string, the str_format which uses to read your tsv file, the type of the string should corresponds exactly with the columns in the tsv file.
                  Defaut parameters, we set these parameters to be some default values, but you can also set it by yourself:
        glm_type: based on the hypothesis, you should define one of the glm types, "group_comparison", "correlation"
        group_label: current group name for this analysis
        freesurfer_home: the environmental variable $FREESURFER_HOME
        surface_file: Specify where to find the data surfaces file in the "CAPS/subject" directory, using specific keywords.
                     For instance, to catch for each subject the cortical thickness, the string used will be :
                     '@subject/@session/t1/freesurfer_cross_sectional/@subject_@session/surf/@hemi.thickness.fwhm@fwhm.fsaverage.mgh'
                     More information is available on the documentation page of the surfstat pipelines. The keywords @subject @ session @hemi @fwhm
                     represents the variable parts.
        path_to_matscript: path to find the matlab script
        full_width_at_half_maximum: fwhm for the surface smoothing, default is 20, integer.
        threshold_uncorrected_pvalue: threshold to display the uncorrected Pvalue, float, default is 0.001.
        threshold_corrected_pvalue: the threshold to display the corrected cluster, default is 0.05, float.
        cluster_threshold: threshold to define a cluster in the process of cluster-wise correction, default is 0.001, float.

    Returns:

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
    clinicasurfstat('%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', %d, '%s', %.3f, '%s', %.3f, '%s', %.3f);
    """ % (input_directory, output_directory, subjects_visits_tsv, design_matrix, contrast, str_format, glm_type, group_label, freesurfer_home, surface_file, 'sizeoffwhm',
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

def json_dict_create(glm_type,
                     design_matrix,
                     str_format,
                     contrast,
                     group_label,
                     full_width_at_half_maximum,
                     threshold_uncorrected_pvalue,
                     threshold_corrected_pvalue,
                     cluster_threshold):
    """
        create a json file containing the glm information
    Args:

        design_matrix: str, the linear model that fits into the GLM, for example '1+group'.
        contrast: string, the contrast matrix for GLM, if the factor you choose is categorized variable, clinica_surfstat will create two contrasts,
                  for example, contrast = 'Label', this will create contrastpos = Label.AD - Label.CN, contrastneg = Label.CN - Label.AD; if the fac-
                  tory that you choose is a continuous factor, clinica_surfstat will just create one contrast, for example, contrast = 'Age', but note,
                  the string name that you choose should be exactly the same with the columns names in your subjects_visits_tsv.
        str_format:string, the str_format which uses to read your tsv file, the typy of the string should corresponds exactly with the columns in the tsv file.
                  Defaut parameters, we set these parameters to be some default values, but you can also set it by yourself:
        glm_type: based on the hypothesis, you should define one of the glm types, "group_comparison", "correlation"
        group_label: current group name for this analysis
        full_width_at_half_maximum: fwhm for the surface smoothing, default is 20, integer.
        threshold_uncorrected_pvalue: threshold to display the uncorrected Pvalue, float, default is 0.001.
        threshold_corrected_pvalue: the threshold to display the corrected cluster, default is 0.05, float.
        cluster_threshold: threshold to define a cluster in the process of cluster-wise correction, default is 0.001, float.

    Returns:

    """
    json_dict = {'AnalysisType': glm_type,
                'DesignMatrix': design_matrix,
                'StringFormatTSV': str_format,
                'Contrast': contrast,
                'GroupLabel': group_label,
                'FWHM': full_width_at_half_maximum,
                'ThresholdUncorrectedPvalue': threshold_uncorrected_pvalue,
                'ThresholdCorrectedPvalue': threshold_corrected_pvalue,
                'ClusterThreshold': cluster_threshold}

    return json_dict


def check_inputs(caps,
                 custom_filename,
                 fwhm,
                 tsv_file):
    """
        Simply checks if the custom strings provided find the correct files. If files are missing, it is easier to
        raise an exception here rather than doing it in the MATLAB script.
    Args:

        caps: The CAPS folder
        custom_filename: The string defining where are the files needed for the surfstat analysis
        fwhm: Full Width at Half Maximum used for the smoothing of the data, needed to catch the file
        tsv_file: tsv file that contains the subject and session lists, along with the other covariates

    Returns:

    """
    import os
    from clinica.utils.stream import cprint
    import pandas as pd

    subjects_visits = pd.io.parsers.read_csv(tsv_file, sep='\t')
    subjects = list(subjects_visits.participant_id)
    sessions = list(subjects_visits.session_id)

    missing_files = []
    for idx in range(len(subjects)):
        fullpath = os.path.join(caps,
                                'subjects',
                                custom_filename.replace('@subject', subjects[idx]).replace('@session', sessions[idx]).replace('@fwhm', str(fwhm)))
        left_hemi = fullpath.replace('@hemi', 'lh')
        right_hemi = fullpath.replace('@hemi', 'rh')

        if not os.path.exists(left_hemi):
            missing_files.append(left_hemi)
        if not os.path.exists(right_hemi):
            missing_files.append(right_hemi)

    if len(missing_files) > 0:
        print(' ** Missing files **')
        for l in missing_files:
            cprint('Not found : ' + l)
        raise Exception(str(len(missing_files)) + ' files not found !')
