from clinica.engine.cworkflow import *

@Visualize("freeview", "-v ${subject_id}/mri/T1.mgz -f ${subject_id}/surf/lh.white:edgecolor=blue ${subject_id}/surf/lh.pial:edgecolor=green ${subject_id}/surf/rh.white:edgecolor=blue ${subject_id}/surf/rh.pial:edgecolor=green", "subject_id")
def recon_all_pipeline(data_dir, output_dir, n_output, recon_all_args='-qcache'):
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
        :param: output_dir: the directory where to put the results of the pipeline
        :param: n_output, scale, the number of output files that you want to contain the results, eg, if you define n_output, then the number of output file should be sub001...sub00(n_output-1)
        :param: recon_all_args, the default value will be set as '-qcache', which will get the result of the fsaverage.
        return: Recon-all workflow
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
    nii_list = []
    for dirpath, dirnames, filenames in os.walk(data_dir):
        subject_list = dirnames
        nii_list = filenames
        break

#    data_dir_out = []
#    for element in subject_list:
#        element = element.split('.')[0]
#        data_dir_out.append(element)
#        
#    data_dir_out = data_dir_out[0: n_output-1]

    wf = pe.Workflow(name='reconall_workflow',base_dir=output_dir)

    datasource = pe.Node(interface = nio.DataGrabber(infields=['subject_id'], outfields=['out_files']), name="datasource")
    datasource.inputs.base_directory = data_dir
    datasource.inputs.template = '%s/%s'
    datasource.inputs.template_args = dict(out_files=[['subject_id', 'nii_list']])
    datasource.inputs.subject_id = subject_list
    datasource.inputs.sort_filelist = True

    try:
        os.makedirs(output_dir)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

    recon_all = pe.MapNode(interface=ReconAll(),name='recon_all', iterfield=['subject_id', 'T1_files'])
    recon_all.inputs.subject_id = subject_list
    recon_all.inputs.subjects_dir = output_dir
    recon_all.inputs.directive = 'all'
    recon_all.inputs.args = recon_all_args;

    wf.connect(datasource,'out_files', recon_all,'T1_files')

    return wf

