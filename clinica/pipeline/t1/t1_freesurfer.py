#!/usr/bin/python#
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 09:04:10 2016

@author: Junhao WEN
"""

from clinica.engine.cworkflow import *

@Visualize("freeview", "-v ${subject_id}/mri/T1.mgz -f ${subject_id}/surf/lh.white:edgecolor=blue ${subject_id}/surf/lh.pial:edgecolor=green ${subject_id}/surf/rh.white:edgecolor=blue ${subject_id}/surf/rh.pial:edgecolor=green", "subject_id")
def recon_all_pipeline(data_dir, output_dir,field_template, template_args, datasink_para = ['orig', 'white'], recon_all_args='-qcache'):
#def recon_all_pipeline(data_dir, output_dir, n_output, recon_all_args='-qcache'):
    
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

        :param: data_dir: the directory where to put the input images, eg, example1.nii, example2.nii, this should be the absolute path
        :param: output_dir: the directory where to put the results of the pipeline, should be absolute path!
        :param: field_template: list, you should define it based on your input data structure       
        :param: template_args: list, you should define it based on your input data structure
        :param: datasink_para: list containing string, optional, the container inside the datasink_folder, for datasinker to store the result that you want, you can define many container to store your result, default is datasink_para = ['orig', 'white']
        :param: recon_all_args, the additional flags for reconAll command line, the default value will be set as '-qcache', which will get the result of the fsaverage.
        return: Recon-all workflow
    """

    import os
    import errno
    import nipype.pipeline.engine as pe
    import nipype.interfaces.io as nio
    from nipype.interfaces.freesurfer.preprocess import ReconAll
    import nipype.interfaces.utility as niu

    try:# just try to check out the ReconAll version for nipype
        if ReconAll.version.fget.func_globals['__version__'].split(".") < ['0', '11', '0']:
            raise RuntimeError('ReconAll version should at least be version of 0.11.0')
    except Exception as e:
        print(str(e))
        exit(1)

    subject_list = []
    # in example_t1_reconall, the data_dir is defined as absolute path, if we should do there, or we should do it here????
    for dirpath, dirnames, filenames in os.walk(data_dir):
        subject_list = dirnames
        break
    output_dir = os.path.expanduser(output_dir)# this is to add ~ before the out_dir, if it doesnt start with ~,just return the path
    try:
        os.makedirs(output_dir)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
    inputnode = pe.MapNode(interface = nio.DataGrabber(infields=['subject_id'], 
                                                     outfields=['out_files']), name="inputnode",
                                                     iterfield = ['subject_id'])   
    recon_all = pe.MapNode(interface=ReconAll(),name='recon_all', iterfield=['subject_id', 'T1_files'])
    outputnode = pe.Node(niu.IdentityInterface(fields=['ReconAll_result']), name='outputnode')
    datasink = pe.Node(nio.DataSink(), name="datasink")
    datasink.inputs.base_directory = output_dir
    datasink.inputs.container = 'datasink_folder'
    
    wf = pe.Workflow(name='reconall_workflow',base_dir=output_dir)
   
    inputnode.inputs.base_directory = data_dir
    inputnode.inputs.template = '*'  
    inputnode.inputs.field_template = dict(out_files=field_template)
    inputnode.inputs.template_args = dict(out_files=[['subject_id', template_args]])

    inputnode.inputs.subject_id = subject_list
    inputnode.inputs.sort_filelist = True

    recon_all.inputs.subject_id = subject_list
    recon_all.inputs.subjects_dir = output_dir
    recon_all.inputs.directive = 'all'
    recon_all.inputs.args = recon_all_args

    wf.connect(inputnode,'out_files', recon_all,'T1_files')
    wf.connect(recon_all, 'subject_id', outputnode, 'ReconAll_result')
    for i in range(len(datasink_para)):
        wf.connect([(recon_all, datasink, [(datasink_para[i], datasink_para[i])])])

    return wf

