def recon_all_pipeline(data_dir,temporary_dir, output_dir):
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

        :param: data_dir: the directory where to put the input images, eg, example1.nii, example2.nii
        :param: temporary_dir: the directory where to put the temporary files in the recon-all pipeline
        :param: output_dir: the directory where to put the results of the pipeline
        :return: Recon-all workflow
    """

    import os
    import errno
    import nipype.pipeline.engine as pe
    import nipype.interfaces.io as nio
    from nipype.interfaces.freesurfer.preprocess import ReconAll

    try:
        if ReconAll.version.fget.func_globals['__version__'].split(".") < ['0', '11', '0']:
            raise RuntimeError('ReconAll version should at least be version of 0.11.0')
    except Exception as e:
        print(str(e))
        exit(1)

    subject_list = []
    for dirpath, dirnames, filenames in os.walk(data_dir):
        subject_list = filenames
        break
    subject_list = filenames

    wf = pe.Workflow(name='reconall_workflow')
    wf.base_dir = temporary_dir

    datasource = pe.Node(interface = nio.DataGrabber(infields=['subject_id'], outfields=['out_files']), name="datasource")
    datasource.inputs.base_directory = data_dir
    datasource.inputs.template = '%s'
    datasource.inputs.subject_id = subject_list
    datasource.inputs.sort_filelist = True

    try:
        os.makedirs(output_dir)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

    data_dir_out = ['sub%03d' % i for i in range(1, 3)]
    recon_all = pe.MapNode(interface=ReconAll(),name='recon_all', iterfield=['subject_id', 'T1_files'])
    recon_all.inputs.subject_id = data_dir_out
    recon_all.inputs.subjects_dir = output_dir
    recon_all.inputs.directive = 'all'

    wf.connect(datasource,'out_files', recon_all,'T1_files')

    return wf

