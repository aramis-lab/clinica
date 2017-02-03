from __future__ import absolute_import
from clinica.pipeline.dwi.dwi_preprocessing import diffusion_preprocessing_t1_based
from clinica.pipeline.dwi.dwi_preprocessing_utils import count_b0s

from os.path import realpath,split,join
import tempfile


data_path = join(split(realpath(__file__))[0], 'external-data/BIDS-example/sub-CLNC01/ses-M00/')

output_directory = tempfile.mkdtemp()
working_directory = tempfile.mkdtemp()

print("Working directory: %s" % working_directory)
print("Output directory: %s" % output_directory)

print("Running...")
dwi_preprocessing = diffusion_preprocessing_t1_based(
    subject_id='ses-CLNC01', session_id='sub-M00', analysis_series_id='default',
    caps_directory=output_directory, num_b0s=count_b0s(join(data_path, 'dwi/sub-CLNC01_ses-M00_dwi.bval')),
    working_directory=working_directory)
dwi_preprocessing.inputs.inputnode.in_t1   = join(data_path, 'anat/sub-CLNC01_ses-M00_T1w.nii.gz')
dwi_preprocessing.inputs.inputnode.in_dwi   = join(data_path, 'dwi/sub-CLNC01_ses-M00_dwi.nii.gz')
dwi_preprocessing.inputs.inputnode.in_bvals  = join(data_path, 'dwi/sub-CLNC01_ses-M00_dwi.bval')
dwi_preprocessing.inputs.inputnode.in_bvecs  = join(data_path, 'dwi/sub-CLNC01_ses-M00_dwi.bvec')
dwi_preprocessing.run()

print("Working directory: %s" % working_directory)
print("Output directory: %s" % output_directory)
