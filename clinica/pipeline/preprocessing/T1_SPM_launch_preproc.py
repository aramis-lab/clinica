from T1_SPM_workflows import T1_SPM_prep_pipeline

import nipype.pipeline.engine as pe
import nipype.interfaces.io as nio
from os import walk


def launch_DataGrabber_example(experiment_dir, output_dir, tissue_map):
    """

    :param experiment_dir:
    :param output_dir:
    :param tissue_map:
    :return:
    """

    data_dir = '/data/FAST_DRIVE2/samper/xnat_download'
    # Subjects
    subjects = []
    for (dirpath, dirnames, filenames) in walk(data_dir):
        subjects = map(lambda x:x[:-4], filenames) #Remove .nii
        break

    #TODO REMOVE
    #subjects = ['005_S_0221','007_S_0068']
    print subjects

    # SelectFiles
    selectfiles = pe.Node(nio.DataGrabber(infields=['subject_id'], outfields=['out_files']), name="selectfiles")
    selectfiles.inputs.base_directory = data_dir
    selectfiles.inputs.template = '%s.nii'
    selectfiles.inputs.subject_id = subjects
    selectfiles.inputs.sort_filelist = False

    T1_SPM_prep_wf = T1_SPM_prep_pipeline(experiment_dir, output_dir, tissue_map)

    preproc_wf = pe.Workflow(name='preproc_wf')
    preproc_wf.base_dir = experiment_dir
    preproc_wf.connect([
        (selectfiles, T1_SPM_prep_wf, [('out_files', 'segmentation_wf.new_segment.channel_files')])
    ])

    preproc_wf.run('MultiProc', plugin_args={'n_procs': 10})


def launch_XNATSource_example(experiment_dir, output_dir, tissue_map):
    """

    :param experiment_dir:
    :param output_dir:
    :param tissue_map:
    :return:
    """

    xnat_input = pe.Node(nio.XNATSource(server='http://134.157.198.180:8080/xnat',
                                        user='george',
                                        pwd='george',
                                        infields=['project', 'subject']),
                         name ='xnat_input')

    xnat_input.inputs.query_template = '/projects/%s/subjects/%s/experiments/*/scans/*/resources/*/files'
    #xnat_input.base_dir = '/data/FAST_DRIVE2/samper/clinica/temp'

    xnat_input.inputs.project = 'test1'
    xnat_input.inputs.subject = ['005_S_0221', '007_S_0068']

    # xnat_input.base_dir = experiment_dir
    # xnat_input.run()

    T1_SPM_prep_wf = T1_SPM_prep_pipeline(experiment_dir, output_dir, tissue_map)

    preproc_wf = pe.Workflow(name='preproc_wf')
    preproc_wf.base_dir = experiment_dir
    preproc_wf.connect([
        (xnat_input, T1_SPM_prep_wf, [('outfiles', 'segmentation_wf.new_segment.channel_files')])
    ])

    preproc_wf.run('MultiProc', plugin_args={'n_procs': 10})


if __name__ == "__main__":

    # Experiment folders
    experiment_dir = '/data/FAST_DRIVE2/samper/clinica/temp'
    output_dir = '/data/FAST_DRIVE2/samper/clinica/temp/output'

    # Experiment specific parameters
    tissue_map = '/aramis/dartagnan2/Software/SPM/spm8_04-2009_updates-r5236-04-02-2013/toolbox/Seg/TPM.nii'

    launch_XNATSource_example(experiment_dir, output_dir, tissue_map)


