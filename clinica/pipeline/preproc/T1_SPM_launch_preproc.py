from T1_SPM_workflows import T1_SPM_prep_pipeline

import nipype.pipeline.engine as pe
import nipype.interfaces.io as nio
from os import walk

# Experiment folders
experiment_dir = '/data/FAST_DRIVE2/samper/clinica'
data_dir = '/data/FAST_DRIVE2/samper/xnat_download'
output_dir = '/data/FAST_DRIVE2/samper/clinica/output'

### Experiment specific parameters
tissue_map = '/aramis/dartagnan2/Software/SPM/spm8_04-2009_updates-r5236-04-02-2013/toolbox/Seg/TPM.nii'

# Subjects
subjects = []
for (dirpath, dirnames, filenames) in walk(data_dir):
    subjects = map(lambda x:x[:-4], filenames) #Remove .nii
    break

#TODO REMOVE
subjects = ['005_S_0221','007_S_0068']
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