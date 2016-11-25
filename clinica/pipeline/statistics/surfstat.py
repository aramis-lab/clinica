#!/usr/bin/python#
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 15:20:40 2016

@author: Junhao WEN
"""

from __future__ import absolute_import
from nipype.interfaces.utility import Function
import nipype.pipeline.engine as pe
import nipype.interfaces.utility as niu
from tempfile import mkdtemp
from clinica.pipeline.statistics.surfstat_utils import absolute_path, get_vars, runmatlab

def clinica_surfstat(input_directory,
                     subjects_visits_tsv,
                     linear_model,
                     contrast,
                     str_format,
                     group_label,
                     analysis_series_id,
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
        Inputs: :param input_directory:  the output folder of recon-all which will contain nested files: ?h.thickness.fwhm**.mgh.
                :param linear_model: string, the linear model that fits into the GLM, for example '1+Lable'.
                :param contrast: string, the contrast matrix for GLM, if the factor you choose is categorized variable, clinica_surfstat will create two contrasts,
                          for example, contrast = 'Label', this will create contrastpos = Label.AD - Label.CN, contrastneg = Label.CN - Label.AD; if the fac-
                          tory that you choose is a continuous factor, clinica_surfstat will just create one contrast, for example, contrast = 'Age', but note,
                          the string name that you choose should be exactly the same with the columns names in your subjects_visits_tsv.
                :param subjects_visits_tsv: string, the path to your tsv file.
                :param analysis_series_id: string, must be the an existed folder(recon-alled output folder.
                :param str_format: string, the str_format which uses to read your tsv file, the typy of the string should corresponds exactly with the columns in the tsv file.
                 Defaut parameters, we set these parameters to be some default values, but you can also set it by yourself:
                :param group_label: Current group name
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
    # Node to fetch the input vars.
    inputnode = pe.Node(name='inputnode',
                        interface=Function(
                            input_names=['input_directory', 'subjects_visits_tsv', 'analysis_series_id', 'group_label'],
                            output_names=['path_to_matscript', 'surfstat_input_dir', 'output_directory'],
                            function=get_vars))
    inputnode.inputs.input_directory = input_directory
    inputnode.inputs.subjects_visits_tsv = subjects_visits_tsv
    inputnode.inputs.analysis_series_id = analysis_series_id
    inputnode.inputs.group_label = group_label

    # Node to wrap the surfstat matlab script.
    surfstat = pe.Node(name='surfstat',
                       interface=Function(input_names=['input_directory', 'output_directory', 'subjects_visits_tsv', 'linear_model',
                                                       'contrast', 'str_format', 'path_to_matscript', 'size_of_fwhm', 'threshold_uncorrected_pvalue',
                                                       'threshold_corrected_pvalue', 'cluster_threshold'],
                                          output_names=['out_images'],
                                          function=runmatlab))
    surfstat.inputs.linear_model = linear_model
    surfstat.inputs.contrast = contrast
    surfstat.inputs.subjects_visits_tsv = subjects_visits_tsv
    surfstat.inputs.str_format = str_format
    surfstat.inputs.size_of_fwhm = size_of_fwhm
    surfstat.inputs.threshold_uncorrected_pvalue = threshold_uncorrected_pvalue
    surfstat.inputs.threshold_corrected_pvalue = threshold_corrected_pvalue
    surfstat.inputs.cluster_threshold = cluster_threshold

    if working_directory is None:
        working_directory = mkdtemp()
    else:
        working_directory = absolute_path(working_directory)

    surfstat_wf = pe.Workflow(name='surfstat_workflow', base_dir=working_directory)

    surfstat_wf.connect(inputnode, 'surfstat_input_dir', surfstat, 'input_directory')
    surfstat_wf.connect(inputnode, 'path_to_matscript', surfstat, 'path_to_matscript')
    surfstat_wf.connect(inputnode, 'output_directory', surfstat, 'output_directory')

    return surfstat_wf
