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
from tempfile import mkdtemp
from clinica.pipeline.t1.t1_freesurfer_utils import create_flags_str, checkfov, absolute_path, write_statistics, log_summary, get_dirs
from clinica.engine.cworkflow import *

@Visualize("freeview", "-v ${subject_id}/mri/T1.mgz -f ${subject_id}/surf/lh.white:edgecolor=blue ${subject_id}/surf/lh.pial:edgecolor=green ${subject_id}/surf/rh.white:edgecolor=blue ${subject_id}/surf/rh.pial:edgecolor=green", "subject_id")
def t1_freesurfer_pipeline(input_dir,
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

        datagrabbernode
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

    # check out ReconAll version
    try:
        if ReconAll.version.fget.func_globals['__version__'].split(".") < ['0', '11', '0']:
            raise RuntimeError('ReconAll version should at least be version of 0.11.0')
    except Exception as e:
        print(str(e))
        exit(1)

    if working_directory is None:
        working_directory = mkdtemp()
    else:
        working_directory = absolute_path(working_directory)

    # Node to get the input vars
    inputnode = pe.Node(name='inputnode',
                          interface=Function(
                          input_names=['output_dir', 'subjects_visits_tsv', 'analysis_series_id'],
                          output_names=['subject_dir', 'subject_id', 'subject_list', 'session_list'],
                          function=get_dirs))
    inputnode.inputs.output_dir = output_dir
    inputnode.inputs.subjects_visits_tsv = subjects_visits_tsv
    inputnode.inputs.analysis_series_id = analysis_series_id

    # Node to grab the BIDS input.
    datagrabbernode = pe.Node(interface=nio.DataGrabber(
                        infields=['subject_list', 'session_list', 'subject_repeat', 'session_repeat'],
                        outfields=['anat_t1']),
                        name="datagrabbernode")  # the best explanation for datagrabber http://nipy.org/nipype/interfaces/generated/nipype.interfaces.io.html#datagrabber
    datagrabbernode.inputs.base_directory = input_dir
    datagrabbernode.inputs.template = '*'
    datagrabbernode.inputs.field_template = dict(anat_t1='%s/%s/anat/%s_%s_T1w.nii.gz')
    # just for test with my home nii files
    # datagrabbernode.inputs.field_template = dict(anat_t1='sub-%s/ses-%s/anat/sub-%s_ses-%s_T1w.nii')
    datagrabbernode.inputs.template_args = dict(anat_t1=[['subject_list', 'session_list', 'subject_repeat',
                                                          'session_repeat']])
    datagrabbernode.inputs.sort_filelist = False

    # MapNode to check out if we need -cw256 for every subject, and -qcache is default for every subject.
    flagnode = pe.MapNode(name='flagnode',
                          iterfield=['t1_list'],
                          interface=Function(
                          input_names=['t1_list', 'recon_all_args'],
                          output_names=['output_flags'],
                          function=checkfov))
    flagnode.inputs.recon_all_args = recon_all_args

    # MapNode to transfer every subject's flag to string.
    create_flags = pe.MapNode(interface=Function(
                              input_names=['input_flags'],
                              output_names=['output_str'],
                              function=create_flags_str),
                              name='create_flags_string',
                              iterfield=['input_flags'])

    # MapNode to implement recon-all.
    recon_all = pe.MapNode(interface=ReconAll(),
                           name='recon_all',
                           iterfield=['subject_id', 'T1_files', 'subjects_dir', 'flags'])
    recon_all.inputs.directive = 'all'

    wf_recon_all = pe.Workflow(name='reconall_workflow')

    wf_recon_all.connect(inputnode, 'subject_list', datagrabbernode, 'subject_list')
    wf_recon_all.connect(inputnode, 'subject_list', datagrabbernode, 'subject_repeat')
    wf_recon_all.connect(inputnode, 'session_list', datagrabbernode, 'session_list')
    wf_recon_all.connect(inputnode, 'session_list', datagrabbernode, 'session_repeat')
    wf_recon_all.connect(inputnode, 'subject_dir', recon_all, 'subjects_dir')
    wf_recon_all.connect(inputnode, 'subject_id', recon_all, 'subject_id')
    wf_recon_all.connect(datagrabbernode, 'anat_t1', recon_all, 'T1_files')
    wf_recon_all.connect(datagrabbernode, 'anat_t1', flagnode, 't1_list')
    wf_recon_all.connect(flagnode, 'output_flags', create_flags, 'input_flags')
    wf_recon_all.connect(create_flags, 'output_str', recon_all, 'flags')

    statisticsnode = pe.Node(name='statisticsnode',
                                interface=Function(
                                 input_names=['subject_dir', 'subject_id', 'analysis_series_id', 'output_dir'],
                                 output_names=['analysis_series_id', 'output_dir'],
                                 function=write_statistics))
    statisticsnode.inputs.analysis_series_id = analysis_series_id
    statisticsnode.inputs.output_dir = output_dir

    lognode = pe.Node(name='lognode',
                      interface=Function(
                          input_names=['subject_list', 'session_list', 'subject_id', 'output_dir', 'analysis_series_id'],
                          output_names=[],
                          function=log_summary))

    wf_recon_all_tsvs = pe.Workflow(name='wf_recon_all_tsvs')

    wf_recon_all_tsvs.connect(statisticsnode, 'analysis_series_id', lognode, 'analysis_series_id')
    wf_recon_all_tsvs.connect(statisticsnode, 'output_dir', lognode, 'output_dir')

    metaflow = pe.Workflow(name='metaflow', base_dir=working_directory)

    metaflow.connect([(wf_recon_all, wf_recon_all_tsvs,[('recon_all.subject_id', 'lognode.subject_id'),
                                                        ('recon_all.subject_id', 'statisticsnode.subject_id'),
                                                        ('recon_all.subject_dir', 'statisticsnode.subject_dir'),
                                                        ('inputnode.subject_list', 'lognode.subject_list'),
                                                        ('inputnode.session_list', 'lognode.session_list'),
                                                        ]),
                       ])

    return metaflow