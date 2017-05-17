#!/usr/bin/python#
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 09:04:10 2016

@author: Junhao Wen
"""
### TODO change the parameters of the name for pipeline and wiki for BIDS

__author__ = "Junhao Wen"
__copyright__ = "Copyright 2016, The Aramis Lab Team"
__credits__ = ["Michael Bacci", "Junhao Wen"]
__license__ = "??"
__version__ = "1.0.0"
__maintainer__ = "Junhao Wen"
__email__ = "junhao.Wen@inria.fr"
__status__ = "Development"


def datagrabber_t1_freesurfer_pipeline(input_dir,
                       output_dir,
                       subjects_visits_tsv=None,
                       working_directory=None,
                       recon_all_args='-qcache'):

    """
        Creates a pipeline for the dataset preparation for "t1_freesurfer_pipeline", this script is based on BIDS dataset,
        if your have non-BIDS dataset which can not use the command line, please adapt them into this pipeline.

        :param: input_dir: str, the path to the directory(default is BIDS) for your dataset.
        :param: output_dir: str, the path to the directory(CAPS) where to put the results of the pipeline.
        :param: subjects_visits_tsv: str, the path pointing to subjects-visit-list tsv file, default behavior is None,
            and this pipeline will run all the subjects in your dataset, if you want to run a sub-group of your dataset,
            please create a separate tsv file for this.
        :param: working_directory: str, path to contain the detail information about your workflow, the default value
            is None, nipype will create a temporary folder to create the pipeline.
        :param: recon_all_args: str, the additional flags for reconAll command line, the default value will be set as
            '-qcache', which will run numerous back-to-back mris_preproc processes for your subjects.

        return: wf_recon_all_with_datagrabber workflow
    """
    import nipype.pipeline.engine as pe
    import errno
    from nipype.interfaces.utility import Function
    from clinica.pipeline.t1.t1_freesurfer_utils import get_dirs_check_reconalled
    from clinica.pipeline.t1.t1_freesurfer_workflows import  t1_freesurfer_pipeline
    from clinica.iotools.utils.data_handling import create_subs_sess_list
    import os
    from clinica.pipeline.t1.t1_freesurfer_utils import bids_datagrabber
    from tempfile import mkdtemp

    if subjects_visits_tsv is None:
        try:
            os.makedirs(output_dir)
        except OSError as exception:
            if exception.errno != errno.EEXIST:  # if the error is not exist error, raise, otherwise, pass
                raise
        tmp_folder = mkdtemp()
        create_subs_sess_list(input_dir, tmp_folder)
        subjects_visits_tsv = os.path.join(tmp_folder, 'subjects_sessions_list.tsv')

    # Node to get the input vars
    inputnode = pe.Node(name='inputnode',
                          interface=Function(
                          input_names=['output_dir', 'subjects_visits_tsv'],
                          output_names=['nouse0', 'nouse1', 'nouse2', 'nouse3', 'subject_dir', 'subject_id', 'subject_list', 'session_list'],
                          function=get_dirs_check_reconalled))
    inputnode.inputs.output_dir = output_dir
    inputnode.inputs.subjects_visits_tsv = subjects_visits_tsv


    # BIDS DataGrabber
    # ===============
    datagrabbernode = pe.Node(name='datagrabbernode',
                              interface=Function(
                                  function=bids_datagrabber,
                                  input_names=['input_dir', 'subjects_list', 'sessions_list'],
                                  output_names=['anat_t1']))
    datagrabbernode.inputs.input_dir = input_dir


    # # Node to grab the BIDS input.
    # datagrabbernode = pe.Node(interface=nio.DataGrabber(
    #                     infields=['subject_list', 'session_list', 'subject_repeat', 'session_repeat'],
    #                     outfields=['anat_t1']),
    #                     name="datagrabbernode")  # the best explanation for datagrabber http://nipy.org/nipype/interfaces/generated/nipype.interfaces.io.html#datagrabber
    # datagrabbernode.inputs.base_directory = input_dir
    # datagrabbernode.inputs.template = '*'
    # datagrabbernode.inputs.field_template = dict(anat_t1='%s/%s/anat/%s_%s_T1w.nii.gz') # This is just to adapt to the output anat_t1
    # # just for test with my home nii files
    # # datagrabbernode.inputs.field_template = dict(anat_t1='sub-%s/ses-%s/anat/sub-%s_ses-%s_T1w.nii')
    # datagrabbernode.inputs.template_args = dict(anat_t1=[['subject_list', 'session_list', 'subject_repeat',
    #                                                       'session_repeat']])
    # datagrabbernode.inputs.sort_filelist = False

    datagrabber_recon_all = t1_freesurfer_pipeline(output_dir,
                           working_directory=working_directory,
                           recon_all_args=recon_all_args)


    wf_recon_all_with_datagrabber = pe.Workflow(name='reconall_workflow', base_dir=working_directory)

    wf_recon_all_with_datagrabber.connect(inputnode, 'subject_list', datagrabbernode, 'subjects_list')
    wf_recon_all_with_datagrabber.connect(inputnode, 'session_list', datagrabbernode, 'sessions_list')
    wf_recon_all_with_datagrabber.connect(inputnode, 'subject_dir', datagrabber_recon_all, 'recon_all.subjects_dir')
    wf_recon_all_with_datagrabber.connect(inputnode, 'subject_id', datagrabber_recon_all, 'recon_all.subject_id')
    wf_recon_all_with_datagrabber.connect(datagrabbernode, 'anat_t1', datagrabber_recon_all, 'recon_all.T1_files')
    wf_recon_all_with_datagrabber.connect(datagrabbernode, 'anat_t1', datagrabber_recon_all, 'flagnode.t1_list')
    wf_recon_all_with_datagrabber.connect(inputnode, 'subject_list', datagrabber_recon_all, 'lognode.subject_list')
    wf_recon_all_with_datagrabber.connect(inputnode, 'session_list', datagrabber_recon_all, 'lognode.session_list')

    return wf_recon_all_with_datagrabber