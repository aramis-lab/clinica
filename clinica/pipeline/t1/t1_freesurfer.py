#!/usr/bin/python#
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 09:04:10 2016

@author: Junhao WEN
"""

import nipype.pipeline.engine as pe
import nipype.interfaces.io as nio
from nipype.interfaces.freesurfer.preprocess import ReconAll
from nipype.interfaces.utility import Function
import nipype.interfaces.utility as niu
from tempfile import mkdtemp
from clinica.pipeline.t1.t1_freesurfer_utils import absolute_path, CAPS_output, write_statistics, checkfov, create_flags_str, get_vars
from clinica.engine.cworkflow import *

@Visualize("freeview", "-v ${subject_id}/mri/T1.mgz -f ${subject_id}/surf/lh.white:edgecolor=blue ${subject_id}/surf/lh.pial:edgecolor=green ${subject_id}/surf/rh.white:edgecolor=blue ${subject_id}/surf/rh.pial:edgecolor=green", "subject_id")
def recon_all_pipeline(input_dir,
                       output_dir,
                       subjects_visits_tsv,
                       analysis_series_id='default',
                       working_directory=None,
                       recon_all_args='-qcache'):

    """
        Creates a pipeline that performs Freesurfer commander, recon-all,
        It takes the input files of MRI T1 images and executes the 31 steps to
        reconstruct the surface of the brain, this progress includes surface-based
        and Volume-based piepeline, which including gray(GM)and white matter(WM)
        segementation, pial and white surface extraction!.

        Inputnode
        ---------
        DataGrabber : FILE
          Mandatory inputs: the input images, should be a string.

        Outputnode
        ----------
        ReconAll
          Optional inputs: T1_files: name of T1 file to process,(a list of items which are an existing file name)
                          args: Additional parameters to the command, (a string)
                          directive: ('all' or 'autorecon1' or 'autorecon2' or 'autorecon2-cp' or 'autorecon2-wm'
                          or 'autorecon2-inflate1' or 'autorecon2-perhemi'
                          or 'autorecon3' or 'localGI' or 'qcache', nipype default value: all)

            For more optional ReconAll inputs and  outputs check:
            http://nipy.org/nipype/interfaces/generated/nipype.interfaces.freesurfer.preprocess.html#reconall

        :param: input_dir: the directory where to put the input images, eg, example1.nii, example2.nii, this should be the absolute path
        :param: output_dir: the directory where to put the results of the pipeline, should be absolute path!
        :param: subjects_visits_tsv: the path pointing to the tsv file
        :param: analysis_series_id: the different rounds that you wanna apply this pipeline for your data, the default value is 'default'
        :param: working_directory: define where to put the infomation of the nipype workflow.
        :param: recon_all_args: the additional flags for reconAll command line, the default value will be set as '-qcache', which will get the result of the fsaverage.

        return: Recon-all workflow
    """


    try:# just try to check out the ReconAll version for nipype
        if ReconAll.version.fget.func_globals['__version__'].split(".") < ['0', '11', '0']:
            raise RuntimeError('ReconAll version should at least be version of 0.11.0')
    except Exception as e:
        print(str(e))
        exit(1)

    subject_dir, subject_id, subject_list, session_list = CAPS_output(output_dir, subjects_visits_tsv, analysis_series_id)

    inputnode = pe.Node(interface=nio.DataGrabber(
                        infields=['subject_id', 'session_id', 'subject_repeat', 'session_repeat'],
                        outfields=['anat_t1']),
                        name="inputnode")  # the best explanation for datagrabber http://nipy.org/nipype/interfaces/generated/nipype.interfaces.io.html#datagrabber
    inputnode.inputs.base_directory = input_dir
    inputnode.inputs.template = '*'
    inputnode.inputs.field_template = dict(anat_t1='%s/%s/anat/%s_%s_T1w.nii.gz')
    # just for test with my home nii files
    # inputnode.inputs.field_template = dict(anat_t1='sub-%s/ses-%s/anat/sub-%s_ses-%s_T1w.nii')
    inputnode.inputs.template_args = dict(anat_t1=[['subject_id', 'session_id', 'subject_repeat', 'session_repeat']]) # the same with dg.inputs.template_args['outfiles']=[['dicomdir','123456-1-1.dcm']]
    inputnode.inputs.subject_id = subject_list
    inputnode.inputs.session_id = session_list
    inputnode.inputs.subject_repeat = subject_list
    inputnode.inputs.session_repeat = session_list
    inputnode.inputs.sort_filelist = False

    flagnode = pe.MapNode(name='flagnode',
                          iterfield=['t1_list'],
                          interface=Function(
                          input_names=['t1_list', 'recon_all_args'],
                          output_names=['output_flags'],
                          function=checkfov))
    flagnode.inputs.recon_all_args = recon_all_args

    create_flags = pe.MapNode(interface=Function(
                              input_names=['input_flags'],
                              output_names=['output_str'],
                              function=create_flags_str),
                              name='create_flags_string',
                              iterfield=['input_flags'])

    recon_all = pe.MapNode(interface=ReconAll(),
                           name='recon_all',
                           iterfield=['subject_id', 'T1_files', 'subjects_dir', 'flags'])
    recon_all.inputs.subject_id = subject_id  # subject_id is the name of every output_subject.
    recon_all.inputs.subjects_dir = subject_dir # subject_dir is the path to contain the different subject folder for the output
    recon_all.inputs.directive = 'all'

    outputnode = pe.Node(niu.IdentityInterface(
                         fields=['ReconAll_result']),
                         name='outputnode')
    if working_directory is None:
        working_directory = mkdtemp()
    else:
        working_directory = absolute_path(working_directory)

    wf_recon_all = pe.Workflow(name='reconall_workflow', base_dir=working_directory)

    wf_recon_all.connect(inputnode, 'anat_t1', recon_all, 'T1_files')
    wf_recon_all.connect(inputnode, 'anat_t1', flagnode, 't1_list')
    wf_recon_all.connect(flagnode, 'output_flags', create_flags, 'input_flags')
    wf_recon_all.connect(create_flags, 'output_str', recon_all, 'flags')
    wf_recon_all.connect(recon_all, 'subject_id', outputnode, 'ReconAll_result')

    return wf_recon_all

def recon_all_statistics_pipeline(output_dir,
                                  subjects_visits_tsv,
                                  analysis_series_id,
                                  working_directory=None):
    """
    Write stats files into tsv files after running recon_all_pipeline() for your dataset.
    :param output_dir:
    :param subjects_visits_tsv:
    :param analysis_series_id: must be the an existed folder(recon-alled output folder)
    :param working_directory:
    :return:
    """

    subject_list, session_list = get_vars(subjects_visits_tsv)

    infosource = pe.Node(niu.IdentityInterface(fields=['subject_id', 'session_id']), name="infosource")
    infosource.inputs.subject_id = subject_list
    infosource.inputs.session_id = session_list

    statisticsnode = pe.MapNode(name='statisticsnode',
                                iterfield=['subject_list', 'session_list'],
                                interface=Function(
                                 input_names=['subject_list', 'session_list', 'analysis_series_id', 'output_dir'],
                                 output_names=[],
                                 function=write_statistics))
    statisticsnode.inputs.analysis_series_id = analysis_series_id
    statisticsnode.inputs.output_dir = output_dir

    if working_directory is None:
        working_directory = mkdtemp()
    else:
        working_directory = absolute_path(working_directory)

    wf_recon_all_statistics = pe.Workflow(name='wf_recon_all_statistics', base_dir=working_directory)
    wf_recon_all_statistics.connect(infosource, 'subject_id', statisticsnode, 'subject_list')
    wf_recon_all_statistics.connect(infosource, 'session_id', statisticsnode, 'session_list')


    return wf_recon_all_statistics