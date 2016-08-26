from t1_spm_workflows import t1_spm_prep_pipeline, segmentation_pipeline

import sys
import argparse
from os import walk

import nipype.pipeline.engine as pe
import nipype.interfaces.io as nio


def launch_segmentation_example(data_dir, experiment_dir, output_dir):
    """
    Example of use of the T1 segmentation workflow with DataGrabber as input.

    :param data_dir: Directory containing the NIFTI files.
    :param experiment_dir: Directory to run the workflow.
    :param output_dir: Directory to save the resulting images of segmentation.
    :return:
    """

    # Retrieving subject list from directory
    subjects = []
    for (dirpath, dirnames, filenames) in walk(data_dir):
        subjects = [x[:-4] for x in filenames if x.endswith('.nii')] #Remove .nii
        break

    # DataGrabber
    selectfiles = pe.Node(nio.DataGrabber(infields=['subject_id'], outfields=['out_files']), name="selectfiles")
    selectfiles.inputs.base_directory = data_dir
    selectfiles.inputs.template = '%s.nii'
    selectfiles.inputs.subject_id = subjects
    selectfiles.inputs.sort_filelist = False

    # Creating T1 preprocessing workflow
    seg_prep_wf = segmentation_pipeline(experiment_dir, output_dir)

    # Creating our own workflow
    seg_wf = pe.Workflow(name='seg_wf')
    seg_wf.base_dir = experiment_dir
    # Connecting the data grabber as input to the preprocessing workflow
    seg_wf.connect([
        (selectfiles, seg_prep_wf, [('out_files', 'new_segment.channel_files')])
    ])

    # Running our workflow
    seg_wf.run('MultiProc', plugin_args={'n_procs': 10})


def launch_DataGrabber_example(data_dir, experiment_dir, output_dir):
    """
    Example of use of the T1 preprocessing workflow with DataGrabber as input.

    :param data_dir: Directory containing the NIFTI files.
    :param experiment_dir: Directory to run the workflow.
    :param output_dir: Directory to save the resulting images of segmentation and registration processes.
    :return:
    """

    # Retrieving subject list from directory
    subjects = []
    for (dirpath, dirnames, filenames) in walk(data_dir):
        subjects = [x[:-4] for x in filenames if x.endswith('.nii')] #Remove .nii
        break

    # DataGrabber
    selectfiles = pe.Node(nio.DataGrabber(infields=['subject_id'], outfields=['out_files']), name="selectfiles")
    selectfiles.inputs.base_directory = data_dir
    selectfiles.inputs.template = '%s.nii'
    selectfiles.inputs.subject_id = subjects
    selectfiles.inputs.sort_filelist = False

    # Creating T1 preprocessing workflow
    t1_spm_prep_wf = t1_spm_prep_pipeline(experiment_dir, output_dir)

    # Creating our own workflow
    preproc_wf = pe.Workflow(name='preproc_wf')
    preproc_wf.base_dir = experiment_dir
    # Connecting the data grabber as input to the preprocessing workflow
    preproc_wf.connect([
        (selectfiles, t1_spm_prep_wf, [('out_files', 'segmentation_wf.new_segment.channel_files')])
    ])

    # Running our workflow
    preproc_wf.run('MultiProc', plugin_args={'n_procs': 10})


def launch_XNATSource_example(server, user, password, experiment_dir, output_dir):
    """
    Example of use of the T1 preprocessing workflow with XNATSource as input.

    :param server: XNAT server address
    :param user: XNAT username
    :param password: XNAT password
    :param experiment_dir: Directory to run the workflow.
    :param output_dir: Directory to save the resulting images of segmentation and registration processes.
    :return:
    """

    # XNATSource recommended initialization is from a configuration file
    # stored under restricted access to keep the database access identifiers safe
    # This is just a dummy example.
    # For more information on XNATSource:
    # http://nipy.org/nipype/interfaces/generated/nipype.interfaces.io.html#xnatsource
    xnat_input = pe.Node(nio.XNATSource(server=server,
                                        user=user,
                                        pwd=password,
                                        infields=['project', 'subject']),
                         name ='xnat_input')

    # A more query can be built using the infields and the query template
    xnat_input.inputs.query_template = '/projects/%s/subjects/%s/experiments/*/scans/*/resources/*/files'

    xnat_input.inputs.project = 'test1'
    xnat_input.inputs.subject = ['005_S_0221', '007_S_0068']

    # Creating T1 preprocessing workflow
    T1_SPM_prep_wf = t1_spm_prep_pipeline(experiment_dir, output_dir)

    # Creating our own workflow
    preproc_wf = pe.Workflow(name='preproc_wf')
    preproc_wf.base_dir = experiment_dir
    # Connecting the XNAT source as input to the preprocessing workflow
    preproc_wf.connect([
        (xnat_input, T1_SPM_prep_wf, [('outfiles', 'segmentation_wf.new_segment.channel_files')])
    ])

    # Running our workflow
    preproc_wf.run('MultiProc', plugin_args={'n_procs': 10})


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Launch SPM T1 pre-processing pipeline')
    parser.add_argument("-d","--data_dir", type=str, help="Path to the directory containing the data to be treated")
    parser.add_argument("-e","--experiment_dir", type=str, help="Path to the directory where to store the experiment data")
    parser.add_argument("-o","--output_dir", type=str, help="Path to the directory where the output should be saved")
    args = parser.parse_args()

    if not args.data_dir or not args.experiment_dir or not args.output_dir:
        print 'Insufficient input arguments'
        parser.print_help()
        sys.exit(1)

    # Experiment folders
    #data_dir = '/data/FAST_DRIVE2/samper/xnat_download'
    #experiment_dir = '/data/FAST_DRIVE2/samper/clinica/run'
    #output_dir = '/data/FAST_DRIVE2/samper/clinica/run/output'


    # To test XNATSource
    # server = 'http://134.157.198.180:8080/xnat'
    # user = 'user'
    # password = 'pwd'

    launch_DataGrabber_example(args.data_dir, args.experiment_dir, args.output_dir)

    sys.exit(0)

    #launch_XNATSource_example(server, user, password, experiment_dir, output_dir)
