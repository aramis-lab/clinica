#!/usr/bin/python

#Import necessary modules from nipype
import os
import errno
from os.path import join as opj
import nipype.pipeline.engine as pe 
import nipype.interfaces.io as nio 
from nipype.interfaces.freesurfer.preprocess import ReconAll

def recon_all_pipeline(data_dir, experiment_dir, output_dir):
    """
    Creates a pipeline that performs Freesurfer commander, recon-all,
    It takes the input files of MRI T1 images and executes the 31 steps to
    reconstruct the surface of the brain, this progress includes surface-based
    and Volume-based piepeline, which including gray(GM)and white matter(WM)
    segementation, pial and white surface extraction!.

    Inputnode
    ---------
    DataGrabber : FILEs
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

    :param: experiment_dir: Directory to run the workflow
    :param: data_dir: Directory to contain the source data to progress
    :param:
    :return: Recon-all workflow
    """
    try:
        if ReconAll.version.fget.func_globals['__version__'].split(".") < ['0','11','0']:
            raise RuntimeError('ReconAll version should at least be version of 0.11.0')
    except Exception as e:
        print(str(e))
        exit(1)
    subject_list = []
    for dirpath, dirnames, filenames in os.walk(data_dir):
        subject_list = filenames
        break
    subject_list = filenames

    #Creat the pipelines
    wf = pe.Workflow(name='reconall_workflow')
    wf.base_dir = opj(data_dir, 'recon-all_working_dir_test')

    #Grab the data
    datasource = pe.Node(interface = nio.DataGrabber(infields=['subject_id'], outfields=['out_files']), name="datasource")
    datasource.inputs.base_directory = data_dir
    datasource.inputs.template = '%s'
    datasource.inputs.subject_id = subject_list
    datasource.inputs.sort_filelist = True

    #make the outputs files
    try:
        os.makedirs(output_dir)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

    #Creat the MapNode
    data_dir_out = ['sub%03d'% i for i in range(1,3)]
    recon_all = pe.MapNode(interface=ReconAll(),name='recon_all', iterfield=['subject_id', 'T1_files'])
    recon_all.inputs.subject_id = data_dir_out
    recon_all.inputs.subjects_dir = output_dir
    recon_all.inputs.directive = 'all'
    recon_all.inputs.args = '-qcache'

    #Connect the Node
    wf.connect(datasource,'out_files', recon_all,'T1_files')
    wf.run("MultiProc", plugin_args={'n_procs':4})
    
    return wf


experiment_dir = '/aramis/home/wen/HAO_lab/Clinica/clinica/examples/data' # the tutorial file to contain all the files
data_dir = opj(experiment_dir, 'Recon-all')
output_dir = opj(data_dir, 'recon-all_test')

print("Data Directory -> %s" % data_dir)
print("Output Directory -> %s" % output_dir)
print("Running...")
recon_all_pipeline(data_dir, experiment_dir, output_dir)

