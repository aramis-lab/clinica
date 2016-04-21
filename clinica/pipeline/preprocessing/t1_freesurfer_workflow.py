#!/usr/bin/python

#Import necessary modules from nipype
import os
from os.path import join as opj
import nipype.pipeline.engine as pe 
import nipype.interfaces.io as nio 
from nipype.interfaces.freesurfer.preprocess import ReconAll

def recon-all_pipeline(data_dir, experiment_dir, output_dir):
"""
    Creates a pipeline that performs Freesurfer commander, recon-all,
    It takes the input files of MRI T1 images and executes the 31 steps to
    reconstruct the surface of the brain, this progress includes surface-based
    and Volume-based piepeline, which including gray(GM)and white matter(WM)
    segementation, pial and white surface extraction!.

    Inputnode
    ---------
    DataGrabber : FILE
      Mandatory input: the input images, should be a string.

    Outputnode
    ----------
    ReconAll
      Optional input: T1_files: name of T1 file to process,(a list of items which are an existing file name)
                      args: Additional parameters to the command, (a string)
                      directive: ('all' or 'autorecon1' or 'autorecon2' or 'autorecon2-cp' or 'autorecon2-wm'
                      or 'autorecon2-inflate1' or 'autorecon2-perhemi'
                      or 'autorecon3' or 'localGI' or 'qcache', nipype default value: all)

        For more optional ReconAll inputs and  outputs check:
        http://nipy.org/nipype/interfaces/generated/nipype.interfaces.freesurfer.preprocess.html#reconall

    :param: experiment_dir: Directory to run the workflow
    :param: data_dir: Directory to contain the source data to progress
    :param:
    :return: Recon-all workflow
    """

    wf = pe.Workflow(name='reconall_workflow')
    wf.base_dir = opj(experiment_dir, 'working_dir')
    wf.base_dir = os.path.abspath('/aramis/home/wen/HAO_lab/Nipype/Data_AD/working_dir')

    datasource = pe.Node(interface = nio.DataGrabber(infields=['subject_id'], outfields=['out_files']), name="datasource")
    datasource.inputs.base_directory = data_dir
    datasource.inputs.template = '%s'
    datasource.inputs.subject_id = subject_list
    datasource.inputs.sort_filelist = True

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    recon_all = pe.MapNode(interface=ReconAll(),name='recon_all', iterfield=['subject_id', 'T1_files'])
    recon_all.inputs.subject_id = data_dir_out
    recon_all.inputs.subjects_dir = output_dir
    recon_all.inputs.directive = 'all'

    wf.connect(datasource,'out_files', recon_all,'T1_files')

    os.rmdir(wf.base_dir)
    return wf
