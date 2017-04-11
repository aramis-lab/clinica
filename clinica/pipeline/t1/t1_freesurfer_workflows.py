#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on Mon Apr 22 09:04:10 2016

@author: Junhao Wen
"""
from clinica.engine.cworkflow import *

__author__ = "Junhao Wen"
__copyright__ = "Copyright 2016, The Aramis Lab Team"
__credits__ = ["Michael Bacci", "Junhao Wen"]
__license__ = "??"
__version__ = "1.0.0"
__maintainer__ = "Junhao Wen"
__email__ = "junhao.Wen@inria.fr"
__status__ = "Development"


@Visualize("freeview", "-v ${subject_id}/mri/T1.mgz -f ${subject_id}/surf/lh.white:edgecolor=blue ${subject_id}/surf/lh.pial:edgecolor=green ${subject_id}/surf/rh.white:edgecolor=blue ${subject_id}/surf/rh.pial:edgecolor=green", "subject_id")
def t1_freesurfer_pipeline(output_dir,
                           working_directory,
                           recon_all_args):
    """
        Creates a pipeline that performs Freesurfer commander, recon-all, It takes the input files of MRI T1 images and
        executes the 31 steps to reconstruct the surface of the brain, this progress includes surface-based and Volume-based
        piepeline, which including gray(GM)and white matter(WM) segementation, pial and white surface extraction!.

        :param: output_dir: str, path to the directory(CAPS) where to put the results of the pipeline.
            default value is 'default'
        :param: working_directory: str, path to contain the detail information about your workflow, the default value
            is None, nipype will create a temporary folder to create the pipeline.
        :param: recon_all_args: str, the additional flags for reconAll command line, the default value will be set as
            '-qcache', which will run numerous back-to-back mris_preproc processes for your subjects.

        return: wf_recon_all

        Note: This workflow is automatically set for BIDS dataset, but if your dataset is not BIDS, please adapt your
        data into this pipeline based on the function "datagrabber_t1_freesurfer_pipeline".
    """
    import nipype.pipeline.engine as pe
    from nipype.interfaces.freesurfer.preprocess import ReconAll
    from nipype.interfaces.utility import Function
    from tempfile import mkdtemp
    from clinica.pipeline.t1.t1_freesurfer_utils import create_flags_str, checkfov, absolute_path, \
        log_summary, write_statistics_per_subject

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

    # MapNode to check out if we need -cw256 for every subject, and -qcache is default for every subject.
    flagnode = pe.MapNode(name='flagnode',
                          iterfield=['t1_list'],
                          interface=Function(
                              input_names=['t1_list', 'recon_all_args'],
                              output_names=['output_flags'],
                              function=checkfov))
    flagnode.inputs.recon_all_args = recon_all_args

    ## TODO: before launch the reconall pipeline, should verify if the number of output_flags equals the number of 'subject_id' from inpurnode, because pybids datagrabber does not return error if the files do not exist


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

    tsvmapnode = pe.MapNode(name='tsvmapnode',
                            iterfield=['subject_id'],
                            interface=Function(
                                input_names=['subject_id', 'output_dir'],
                                output_names=[],
                                function=write_statistics_per_subject))
    tsvmapnode.inputs.output_dir = output_dir

    lognode = pe.Node(name='lognode',
                      interface=Function(
                          input_names=['subject_list', 'session_list', 'subject_id', 'output_dir'],
                          output_names=[],
                          function=log_summary))
    lognode.inputs.output_dir = output_dir

    wf_recon_all = pe.Workflow(name='reconall_workflow', base_dir=working_directory)

    wf_recon_all.connect(flagnode, 'output_flags', create_flags, 'input_flags')
    wf_recon_all.connect(create_flags, 'output_str', recon_all, 'flags')
    wf_recon_all.connect(recon_all, 'subject_id', tsvmapnode, 'subject_id')
    wf_recon_all.connect(recon_all, 'subject_id', lognode, 'subject_id')

    return wf_recon_all