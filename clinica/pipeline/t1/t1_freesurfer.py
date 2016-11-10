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

    def absolute_path(arg):
        """Transfer any path to absolute path"""
        if arg[:1] == '~':
            return os.path.expanduser(arg)
        elif arg[:1] == '.':
            return os.getcwd()
        else:
            return os.path.join(os.getcwd(), arg)

    def CAPS_output(output_dir):
        """Define and create the CAPS output """
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

    subject_dir, subject_id, subject_list, session_list = CAPS_output(output_dir)

    def checkFOV(t1_list, recon_all_args):
        """Verifying size of inputs and FOV of each T1 image"""

        import sys
        import nibabel as nib
        from nipype.utils.filemanip import filename_to_list

        output_flags = []
        num_t1 = len(t1_list)
        t1_list = filename_to_list(t1_list)

        if num_t1 == 0:
            print("ERROR: No T1's Given")
            sys.exit(-1)

        shape = nib.load(t1_list[0]).shape
        for t1 in t1_list:
            f = nib.load(t1)
            voxel_size = f.header.get_zooms()
            t1_size = f.header.get_data_shape()
            # not sure if we should constrain all the T1 file should have the same size
            # if t1_size != shape:
            #     print("ERROR: T1s not the same size. Cannot process {0} and {1} "
            #           "together".format(t1_list[0], t1))
            #     sys.exit(-1)
            if voxel_size[0] * t1_size[0] > 256 or voxel_size[1] * t1_size[1] or voxel_size[2] * t1_size[2]:
                print("Setting MRI Convert to crop images to 256 FOV")
                optional_flag = '-cw256'
            else:
                print("No need to add -cw256 flag")
                optional_flag = ''
            flag = "{0} ".format(recon_all_args) + optional_flag
            output_flags.append(flag)

        return output_flags

    def create_flags_str(input_flags):
        """
        Create a commandline string from a list of input flags
        """
        output_str = ""
        for flag in input_flags:
            output_str += "{0} ".format(flag)
        output_str.strip() #stripped from the beginning and the end of the string (default whitespace characters).

        return output_str

    inputnode = pe.Node(interface=nio.DataGrabber(
                        infields=['subject_id', 'session', 'subject_repeat', 'session_repeat'],
                        outfields=['anat_t1']),
                        name="inputnode")  # the best explanation for datagrabber http://nipy.org/nipype/interfaces/generated/nipype.interfaces.io.html#datagrabber

    internode = pe.MapNode(name='internode',
                           iterfield=['t1_list'],
                           interface=Function(
                           input_names=['t1_list', 'recon_all_args'],
                           output_names=['output_flags'], function=checkFOV))
    internode.inputs.recon_all_args = recon_all_args

    create_flags = pe.MapNode(interface=Function(
                              input_names=['input_flags'],
                              output_names=['output_str'],
                              function=create_flags_str),
                              name='create_flags_string',
                              iterfield=['input_flags'])

    recon_all = pe.MapNode(interface=ReconAll(),
                           name='recon_all',
                           iterfield=['subject_id', 'T1_files', 'subjects_dir', 'flags'])

    outputnode = pe.Node(niu.IdentityInterface(
                         fields=['ReconAll_result']),
                         name='outputnode')

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

    wf.connect(inputnode, 'anat_t1', recon_all, 'T1_files')
    wf.connect(inputnode, 'anat_t1', internode, 't1_list')
    wf.connect(internode, 'output_flags', create_flags, 'input_flags')
    wf.connect(create_flags, 'output_str', recon_all, 'flags')
    wf.connect(recon_all, 'subject_id', outputnode, 'ReconAll_result')

    return wf
