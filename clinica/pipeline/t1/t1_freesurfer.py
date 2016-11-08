#!/usr/bin/python#
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 09:04:10 2016

@author: Junhao WEN
"""

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

    import os, csv
    import errno
    import nipype.pipeline.engine as pe
    import nipype.interfaces.io as nio
    from nipype.interfaces.freesurfer.preprocess import ReconAll
    from nipype.interfaces.utility import Function
    import nipype.interfaces.utility as niu
    from tempfile import mkdtemp

    try:# just try to check out the ReconAll version for nipype
        if ReconAll.version.fget.func_globals['__version__'].split(".") < ['0', '11', '0']:
            raise RuntimeError('ReconAll version should at least be version of 0.11.0')
    except Exception as e:
        print(str(e))
        exit(1)

    #transfer any path to be absolute path.
    def absolute_path(arg):
        if arg[:1] == '~':
            return os.path.expanduser(arg)
        elif arg[:1] == '.':
            return os.getcwd()
        else:
            return os.path.join(os.getcwd(), arg)

    # new version for BIDS
    def BIDS_output(output_dir):
        # last_dir = os.path.basename(input_dir)
        # dataset_name = last_dir.split('_')[0]

        subject_list = []
        session_list = []
        subject_id=[]
        with open(subjects_visits_tsv, 'rb') as tsvin:
            tsv_reader = csv.reader(tsvin, delimiter='\t')

            for row in tsv_reader:
                subject_list.append(row[0])
                session_list.append(row[1])
                subject_id.append(row[0] + '_' + row[1])

        output_path = os.path.expanduser(output_dir)  # this is to add ~ before the out_dir, if it doesnt start with ~,just return the path
        output_base = 'analysis-series-' + analysis_series_id + '/subjects'
        output_dir = output_path + '/' + output_base
        try:
            os.makedirs(output_dir)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise

        subject_dir = []
        num_subject = len(subject_list)
        for i in range(num_subject):
            subject = output_dir + '/' + 'sub-' + subject_list[i] + '/' + 'ses-' + session_list[i] + '/' + 't1' + '/' + 'freesurfer-cross-sectional'
            try:
                os.makedirs(subject)
            except OSError as exception:
                if exception.errno != errno.EEXIST:
                    raise
            subject_dir.append(subject)

        return subject_dir, subject_id, subject_list, session_list

    subject_dir, subject_id, subject_list, session_list = BIDS_output(output_dir)

    #
    # def checkoutFOV(t1_list):
    #     import nibabel as nib
    #     output_flag = []
    #     for i in t1_list:
    #         f = nib.load(i)
    #         voxel_size = f.header.get_zooms()
    #         data_shape = f.header.get_data_shape()
    #         if voxel_size[0] * data_shape[0] > 256 or voxel_size[1] * data_shape[1] or voxel_size[2] * data_shape[2]:
    #             output_flag.append('cw256')
    #         else:
    #             output_flag.append('all')
    #     print(output_flag)
    #     return output_dir
    #
    # def checkFOV(T1_files):
    #     """Verifying size of inputs"""
    #
    #     # modified from nipype.workflows.smri.freesurfer.autorecon1.checkT1s:
    #     # https://github.com/nipy/nipype/blob/master/nipype/workflows/smri/freesurfer/autorecon1.py
    #
    #     import sys
    #     import nibabel as nib
    #     from nipype.utils.filemanip import filename_to_list
    #
    #     output_flags = []
    #
    #     T1_files = filename_to_list(T1_files)
    #     if len(T1_files) == 0:
    #         print("ERROR: No T1's Given")
    #         sys.exit(-1)
    #
    #     shape = nib.load(T1_files[0]).shape
    #     for t1 in T1_files[1:]:
    #         if nib.load(t1).shape != shape:
    #             print("ERROR: T1s not the same size. Cannot process {0} and {1} "
    #                   "together".format(T1_files[0], t1))
    #             sys.exit(-1)
    #
    #     # check if cw256 is set to crop the images if size is larger than 256
    #     if any(dim > 256 for dim in shape):
    #         print("Setting MRI Convert to crop images to 256 FOV")
    #         output_flags.append('-cw256')
    #     else:
    #         print("No need to add -cw256 flag")
    #         output_flags.append(None)
    #
    #     return output_flags
    #
    inputnode = pe.Node(interface=nio.DataGrabber(infields=['subject_id', 'session', 'subject_repeat', 'session_repeat'],
                                                outfields=['anat_t1']), name="inputnode")  # the best explanation for datagrabber http://nipy.org/nipype/interfaces/generated/nipype.interfaces.io.html#datagrabber
    #
    # internode = pe.Node(name='internode',
    #                     interface=Function(
    #                        input_names=['t1_list'],
    #                        output_names=['output_flag'], function=checkFOV))

    recon_all = pe.MapNode(interface=ReconAll(), name='recon_all', iterfield=['subject_id', 'T1_files', 'subjects_dir'])
    # create a special identity interface for outputing the subject_id
    outputnode = pe.Node(niu.IdentityInterface(fields=['ReconAll_result']), name='outputnode')
    #datasink = pe.Node(nio.DataSink(), name="datasink")
    #datasink.inputs.base_directory = output_dir
    #datasink.inputs.container = 'datasink_folder'

    # define the base_dir to get the reconall_workflow info, if not set, default=None, which results in the use of mkdtemp
    # wf = pe.Workflow(name='reconall_workflow', base_dir=output_dir)
    # if the user want to keep this workflow infor, this can be set as an optional prarm, but if workflow is stopped manually, not finished, the workflow info seems to be in the current folder, should test it.
    if working_directory is None:
        working_directory = mkdtemp()
    else:
        working_directory = absolute_path(working_directory)

    wf = pe.Workflow(name='reconall_workflow', base_dir=working_directory)



    inputnode.inputs.base_directory = input_dir
    inputnode.inputs.template = '*'
    inputnode.inputs.field_template = dict(anat_t1='sub-%s/ses-%s/anat/sub-%s_ses-%s_T1w.nii.gz')
    # just for test with my home nii files
    # inputnode.inputs.field_template = dict(anat_t1='sub-%s/ses-%s/anat/sub-%s_ses-%s_T1w.nii')

    inputnode.inputs.template_args = dict(anat_t1=[['subject_id', 'session', 'subject_repeat', 'session_repeat']]) # the same with dg.inputs.template_args['outfiles']=[['dicomdir','123456-1-1.dcm']]
    inputnode.inputs.subject_id = subject_list
    inputnode.inputs.session = session_list
    inputnode.inputs.subject_repeat = subject_list
    inputnode.inputs.session_repeat = session_list
    inputnode.inputs.sort_filelist = False


    recon_all.inputs.subject_id = subject_id  # subject_id is the name of every output_subject.
    recon_all.inputs.subjects_dir = subject_dir # subject_dir is the path to contain the different subject folder for the output
    recon_all.inputs.directive = 'all'
    recon_all.inputs.args = recon_all_args
    # recon_all.inputs.flags = internode.outputs.output_flag
    # print("output flag are: %s" % internode.outputs.output_flag)

    wf.connect(inputnode, 'anat_t1', recon_all, 'T1_files')
    # wf.connect(inputnode, 'anat_t1', internode, 't1_list')
    wf.connect(recon_all, 'subject_id', outputnode, 'ReconAll_result')
    #for i in range(len(datasink_para)):
        #wf.connect([(recon_all, datasink, [(datasink_para[i], datasink_para[i])])])

    return wf
