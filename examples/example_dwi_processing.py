#!/usr/bin/python
"""This module launches a tractography pipeline."""

from __future__ import absolute_import
from clinica.pipeline.dwi.dwi_processing import dwi_processing_pipeline

from os.path import join

#path_to_data = join(split(realpath(__file__))[0], 'data/subject_example')

data_path = '/Users/alexandre.routier/Data/subject_example'
output_directory = "/tmp/subject_example"

dwi_processing = dwi_processing_pipeline(
    participant_id='CLNC01', session_id='M00',
    caps_directory=output_directory, tractography_nb_of_tracks="100K" , nthreads=2)
dwi_processing.inputs.inputnode.in_dwi_nii                  = join(data_path, 'dwi_preprocessing/preprocessing/out_file/vol0000_maths_thresh_merged.nii.gz')
dwi_processing.inputs.inputnode.in_bvecs                    = join(data_path, 'dwi_preprocessing/preprocessing/out_bvecs/bvecs_rotated.bvec')
dwi_processing.inputs.inputnode.in_bvals                    = join(data_path, 'dwi_preprocessing/preprocessing/out_bval/bvals')
dwi_processing.inputs.inputnode.in_b0_mask                  = join(data_path, 'dwi_preprocessing/preprocessing/b0_mask/vol0000_warp_maths_thresh_merged_roi_brain_mask.nii.gz')
dwi_processing.inputs.inputnode.in_white_matter_binary_mask = join(data_path, 'fsl_t1_segmentation/out_tissue_class_files/fast_seg_2.nii.gz')



print("Running tractography pipeline...")
dwi_processing.run()
print("...Over!")
