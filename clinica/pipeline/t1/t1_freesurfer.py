#!/usr/bin/python#
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 09:04:10 2016

@author: Junhao WEN
"""

from clinica.engine.cworkflow import *

@Visualize("freeview", "-v ${subject_id}/mri/T1.mgz -f ${subject_id}/surf/lh.white:edgecolor=blue ${subject_id}/surf/lh.pial:edgecolor=green ${subject_id}/surf/rh.white:edgecolor=blue ${subject_id}/surf/rh.pial:edgecolor=green", "subject_id")
def recon_all_pipeline(data_dir, output_dir, tsv_file, recon_all_args='-qcache'):

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
        :param: tsv_file: the path pointing to the tsv file
        :param: recon_all_args, the additional flags for reconAll command line, the default value will be set as '-qcache', which will get the result of the fsaverage.
        return: Recon-all workflow
    """

    import os, csv
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

# new version for BIDS
    def BIDS_input(data_dir, output_dir):
        base_dir = os.path.basename(data_dir)
        dataset_name = base_dir.split('_')[0]

        subject_list = []
        session_list = []
        subject_session=[]
        with open(tsv_file, 'rb') as tsvin:
            tsv_reader = csv.reader(tsvin, delimiter='\t')

            for row in tsv_reader:
                subject_list.append(row[0])
                session_list.append(row[1])
                subject_session.append(row[0] + '_' + row[1])

        output_path = os.path.expanduser(output_dir)  # this is to add ~ before the out_dir, if it doesnt start with ~,just return the path
        output_base = dataset_name + '_CAPS/analysis-series-default/subjects'
        output_dir = output_path + '/' + output_base
        try:
            os.makedirs(output_dir)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise
        return subject_list, session_list, output_dir, subject_session

    subject_list, session_list, output_dir, subject_session = BIDS_input(data_dir, output_dir)

    def BIDS_output(output_dir):
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
        return subject_dir

    inputnode = pe.Node(interface=nio.DataGrabber(infields=['subject_id', 'session', 'subject_repeat', 'session_repeat'],
                                                outfields=['anat_t1']), name="inputnode")  # the best explanation for datagrabber http://nipy.org/nipype/interfaces/generated/nipype.interfaces.io.html#datagrabber

    recon_all = pe.MapNode(interface=ReconAll(), name='recon_all', iterfield=['subject_id', 'T1_files', 'subjects_dir'])
    outputnode = pe.Node(niu.IdentityInterface(fields=['ReconAll_result']), name='outputnode')
    #datasink = pe.Node(nio.DataSink(), name="datasink")
    #datasink.inputs.base_directory = output_dir
    #datasink.inputs.container = 'datasink_folder'

    # define the base_dir to get the reconall_workflow info, if not set, default=None, which results in the use of mkdtemp
    # wf = pe.Workflow(name='reconall_workflow', base_dir=output_dir)
    wf = pe.Workflow(name='reconall_workflow')



    inputnode.inputs.base_directory = data_dir
    inputnode.inputs.template = '*'
    inputnode.inputs.field_template = dict(anat_t1='sub-%s/ses-%s/anat/sub-%s_ses-%s_T1w.nii.gz')
    inputnode.inputs.template_args = dict(anat_t1=[['subject_id', 'session', 'subject_repeat', 'session_repeat']]) # the same with dg.inputs.template_args['outfiles']=[['dicomdir','123456-1-1.dcm']]
    inputnode.inputs.subject_id = subject_list
    inputnode.inputs.session = session_list
    inputnode.inputs.subject_repeat = subject_list
    inputnode.inputs.session_repeat = session_list
    inputnode.inputs.sort_filelist = False


    recon_all.inputs.subject_id = subject_session  # subject_id is the subject folder name of the output of reconall.
    recon_all.inputs.subjects_dir = BIDS_output(output_dir) # subject_dir is the path to contain the different subject folder for the output
    recon_all.inputs.directive = 'all'
    recon_all.inputs.args = recon_all_args

    wf.connect(inputnode, 'anat_t1', recon_all, 'T1_files')
    wf.connect(recon_all, 'subject_id', outputnode, 'ReconAll_result')
    #for i in range(len(datasink_para)):
        #wf.connect([(recon_all, datasink, [(datasink_para[i], datasink_para[i])])])

    return wf
