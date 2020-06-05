# coding: utf8

"""
This file contains a set of functional tests designed to check the correct execution of the pipeline and the
different functions available in Clinica
"""

__author__ = "Arnaud Marcoux"
__copyright__ = "Copyright 2016-2019 The Aramis Lab Team"
__credits__ = ["Arnaud Marcoux"]
__license__ = "See LICENSE.txt file"
__version__ = "0.2.0"
__maintainer__ = "Arnaud Marcoux, Mauricio Diaz"
__email__ = "arnaud.marcoux@inria.fr, mauricio.diaz@inria.fr"
__status__ = "Development"


import warnings
from os import pardir

from testing_tools import *

# Determine location for working_directory
warnings.filterwarnings("ignore")


def test_run_T1FreeSurferCrossSectional(cmdopt):
    # Data for this functional test comes from https://openneuro.org/datasets/ds000204
    from os.path import dirname, join, abspath
    from clinica.pipelines.t1_freesurfer.t1_freesurfer_pipeline import T1FreeSurfer

    working_dir = cmdopt
    root = dirname(abspath(join(abspath(__file__), pardir)))
    root = join(root, 'data', 'T1FreeSurfer')

    # Remove potential residual of previous tests
    clean_folder(join(root, 'out', 'caps'), recreate=False)
    clean_folder(join(working_dir, 'T1FreeSurfer'))

    parameters = {
        'recon_all_args': '-qcache',
    }

    pipeline = T1FreeSurfer(
        bids_directory=join(root, 'in', 'bids'),
        caps_directory=join(root, 'out', 'caps'),
        tsv_file=join(root, 'in', 'subjects.tsv'),
        parameters=parameters,
        base_dir=join(working_dir, 'T1FreeSurfer'),
    )
    pipeline.base_dir = join(working_dir, 'T1FreeSurfer')
    pipeline.run(bypass_check=True)

    # We only check that folders are the same meaning that FreeSurfer finished without error
    # surf/ folder is ignored because it contains sym links that makes hard to check with ref data
    # (sym links of ref data are ignored after rsync on CI machines)
    def path_to_caps_fs(part_id, sess_id):
        import os
        output_folder = os.path.join('caps', 'subjects', part_id, sess_id, 't1', 'freesurfer_cross_sectional')
        return output_folder

    compare_folders(join(root, 'out'), join(root, 'ref'),
                    join(path_to_caps_fs('sub-01', 'ses-2011'), 'regional_measures'))
    compare_folders(join(root, 'out'), join(root, 'ref'),
                    join(path_to_caps_fs('sub-01', 'ses-2011'), 'sub-01_ses-2011', 'label'))
    compare_folders(join(root, 'out'), join(root, 'ref'),
                    join(path_to_caps_fs('sub-01', 'ses-2011'), 'sub-01_ses-2011', 'mri'))
    compare_folders(join(root, 'out'), join(root, 'ref'),
                    join(path_to_caps_fs('sub-01', 'ses-2011'), 'sub-01_ses-2011', 'stats'))

    clean_folder(join(root, 'out', 'caps'), recreate=False)
    clean_folder(join(working_dir, 'T1FreeSurfer'), recreate=False)


def test_run_T1VolumeTissueSegmentation(cmdopt):
    import os
    from os.path import dirname, join, abspath
    from clinica.pipelines.t1_volume_tissue_segmentation.t1_volume_tissue_segmentation_pipeline import T1VolumeTissueSegmentation

    working_dir = cmdopt
    root = dirname(abspath(join(abspath(__file__), pardir)))
    root = join(root, 'data', 'T1VolumeTissueSegmentation')
    clean_folder(join(working_dir, 'T1VolumeTissueSegmentation'))
    clean_folder(join(root, 'out', 'caps'))

    pipeline = T1VolumeTissueSegmentation(
        bids_directory=join(root, 'in', 'bids'),
        caps_directory=join(root, 'out', 'caps'),
        tsv_file=join(root, 'in', 'subjects.tsv'),
        base_dir=join(working_dir, 'T1VolumeTissueSegmentation'),
    )
    pipeline.build()
    pipeline.run(bypass_check=True)

    out_file = join(root, 'out/caps/subjects/sub-ADNI011S4105/ses-M00/t1/spm/segmentation/dartel_input/'
                    + 'sub-ADNI011S4105_ses-M00_T1w_segm-graymatter_dartelinput.nii.gz')
    if not os.path.exists(out_file):
        raise IOError('Pipeline did not produce file: ' + out_file + '. Consider rerunning test_run_T1VolumeTissueSegmentation')

    ref_file = join(root, 'ref/caps/subjects/sub-ADNI011S4105/ses-M00/t1/spm/segmentation/dartel_input/'
                    + 'sub-ADNI011S4105_ses-M00_T1w_segm-graymatter_dartelinput.nii.gz')

    assert likeliness_measure(out_file, ref_file, (1e-1, 0.02), (0.4, 0.01))

    clean_folder(join(root, 'out', 'caps'), recreate=False)
    clean_folder(join(working_dir, 'T1VolumeTissueSegmentation'), recreate=False)


def test_run_T1VolumeCreateDartel(cmdopt):
    from os.path import dirname, join, abspath
    import shutil
    from clinica.pipelines.t1_volume_create_dartel.t1_volume_create_dartel_pipeline import T1VolumeCreateDartel

    working_dir = cmdopt
    root = dirname(abspath(join(abspath(__file__), pardir)))
    root = join(root, 'data', 'T1VolumeCreateDartel')

    # Remove potential residual of previous UT
    clean_folder(join(working_dir, 'T1VolumeCreateDartel'))
    clean_folder(join(root, 'out', 'caps'), recreate=False)
    shutil.copytree(join(root, 'in', 'caps'), join(root, 'out', 'caps'))

    parameters = {
        'group_id': 'UnitTest'
    }
    # Instantiate pipeline
    pipeline = T1VolumeCreateDartel(
        bids_directory=join(root, 'in', 'bids'),
        caps_directory=join(root, 'out', 'caps'),
        tsv_file=join(root, 'in', 'subjects.tsv'),
        base_dir=join(working_dir, 'T1VolumeCreateDartel'),
        parameters=parameters
    )
    pipeline.build()
    pipeline.run(plugin='MultiProc', plugin_args={'n_procs': 4}, bypass_check=True)

    # Check output vs ref
    out_data_template = join(root, 'out/caps/groups/group-UnitTest/t1/group-UnitTest_template.nii.gz')
    ref_data_template = join(root, 'ref/group-UnitTest_template.nii.gz')
    assert likeliness_measure(out_data_template, ref_data_template, (1e-3, 0.1), (1e-2, 0.1))

    subjects = ['sub-ADNI011S4105', 'sub-ADNI023S4020', 'sub-ADNI035S4082', 'sub-ADNI128S4832']
    out_data_forward_def = [join(root, 'out', 'caps', 'subjects', sub, 'ses-M00', 't1', 'spm', 'dartel',
                                 'group-UnitTest', sub +
                                 '_ses-M00_T1w_target-UnitTest_transformation-forward_deformation.nii.gz')
                            for sub in subjects]
    ref_data_forward_def = [join(root, 'ref', sub
                                 + '_ses-M00_T1w_target-UnitTest_transformation-forward_deformation.nii.gz')
                            for sub in subjects]

    for i in range(len(out_data_forward_def)):
        assert likeliness_measure(out_data_forward_def[i], ref_data_forward_def[i], (1e-3, 0.25), (1e-2, 0.1))

    # Remove data in out folder
    clean_folder(join(root, 'out', 'caps'), recreate=False)
    clean_folder(join(working_dir, 'T1VolumeCreateDartel'), recreate=False)


def test_run_T1VolumeDartel2MNI(cmdopt):
    from os.path import dirname, join, abspath
    import shutil
    from clinica.pipelines.t1_volume_dartel2mni.t1_volume_dartel2mni_pipeline import T1VolumeDartel2MNI

    working_dir = cmdopt
    root = dirname(abspath(join(abspath(__file__), pardir)))
    root = join(root, 'data', 'T1VolumeDartel2MNI')

    # Remove potential residual of previous UT
    clean_folder(join(root, 'out', 'caps'), recreate=False)
    clean_folder(join(working_dir, 'T1VolumeDartel2MNI'))

    # Copy necessary data from in to out
    shutil.copytree(join(root, 'in', 'caps'), join(root, 'out', 'caps'))

    parameters = {
        'group_id': 'UnitTest'
    }
    # Instantiate pipeline and run()
    pipeline = T1VolumeDartel2MNI(
        bids_directory=join(root, 'in', 'bids'),
        caps_directory=join(root, 'out', 'caps'),
        tsv_file=join(root, 'in', 'subjects.tsv'),
        base_dir=join(working_dir, 'T1VolumeDartel2MNI'),
        parameters=parameters
    )
    pipeline.build()
    pipeline.run(plugin='MultiProc', plugin_args={'n_procs': 4}, bypass_check=True)

    # Check output vs ref
    subjects = ['sub-ADNI011S4105', 'sub-ADNI023S4020', 'sub-ADNI035S4082', 'sub-ADNI128S4832']
    out_data_GM_MNI = [join(root, 'out', 'caps', 'subjects', sub, 'ses-M00', 't1', 'spm', 'dartel', 'group-UnitTest',
                            sub + '_ses-M00_T1w_segm-graymatter_space-Ixi549Space_modulated-on_fwhm-8mm_probability.nii.gz')
                       for sub in subjects]
    ref_data_GM_MNI = [join(root, 'ref', sub + '_ses-M00_T1w_segm-graymatter_space-Ixi549Space_modulated-on_fwhm-8mm_probability.nii.gz')
                       for sub in subjects]
    for i in range(len(out_data_GM_MNI)):
        assert likeliness_measure(out_data_GM_MNI[i], ref_data_GM_MNI[i],
                                  (1e-4, 0.15), (1, 0.02))

    # Remove data in out folder
    clean_folder(join(root, 'out', 'caps'), recreate=False)
    clean_folder(join(working_dir, 'T1VolumeDartel2MNI'), recreate=False)


def test_run_T1VolumeRegisterDartel(cmdopt):
    from os.path import dirname, join, abspath
    import shutil
    from clinica.pipelines.t1_volume_register_dartel.t1_volume_register_dartel_pipeline import T1VolumeRegisterDartel

    working_dir = cmdopt
    root = dirname(abspath(join(abspath(__file__), pardir)))
    root = join(root, 'data', 'T1VolumeExistingDartel')
    clean_folder(join(root, 'out', 'caps'), recreate=False)
    clean_folder(join(working_dir, 'T1VolumeExistingDartel'))

    # Copy necessary data to run pipeline
    shutil.copytree(join(root, 'in', 'caps'), join(root, 'out', 'caps'))

    # Instantiate and run pipeline
    parameters = {
        'group_id': 'UnitTest'
    }
    pipeline = T1VolumeRegisterDartel(
        bids_directory=join(root, 'in', 'bids'),
        caps_directory=join(root, 'out', 'caps'),
        tsv_file=join(root, 'in', 'subjects.tsv'),
        base_dir=join(working_dir, 'T1VolumeExistingDartel'),
        parameters=parameters
    )
    pipeline.build()
    pipeline.run(plugin='MultiProc', plugin_args={'n_procs': 4}, bypass_check=True)

    # Check output vs ref
    subjects = ['sub-ADNI011S4105', 'sub-ADNI023S4020', 'sub-ADNI035S4082', 'sub-ADNI128S4832']
    out_data_forward_def = [join(root, 'out', 'caps', 'subjects', sub, 'ses-M00', 't1', 'spm', 'dartel',
                                 'group-UnitTest',
                                 sub + '_ses-M00_T1w_target-UnitTest_transformation-forward_deformation.nii.gz')
                            for sub in subjects]
    ref_data_forward_def = [join(root, 'ref',
                                 sub + '_ses-M00_T1w_target-UnitTest_transformation-forward_deformation.nii.gz')
                            for sub in subjects]

    for i in range(len(out_data_forward_def)):
        assert likeliness_measure(out_data_forward_def[i], ref_data_forward_def[i], (1e-3, 0.25), (1e-2, 0.1))

    # Remove data in out folder
    clean_folder(join(root, 'out', 'caps'), recreate=False)
    clean_folder(join(working_dir, 'T1VolumeExistingDartel'), recreate=False)


def test_run_T1VolumeParcellation(cmdopt):
    from os.path import dirname, join, abspath
    import shutil
    import pandas as pds
    import numpy as np
    from clinica.pipelines.t1_volume_parcellation.t1_volume_parcellation_pipeline import T1VolumeParcellation

    working_dir = cmdopt
    root = dirname(abspath(join(abspath(__file__), pardir)))
    root = join(root, 'data', 'T1VolumeParcellation')
    clean_folder(join(root, 'out', 'caps'), recreate=False)
    clean_folder(join(working_dir, 'T1VolumeParcellation'))

    # Copy data for use of pipeline
    shutil.copytree(join(root, 'in', 'caps'), join(root, 'out', 'caps'))

    # Instantiate pipeline
    parameters = {
        'group_id': 'UnitTest'
    }
    pipeline = T1VolumeParcellation(
        caps_directory=join(root, 'in', 'caps'),
        tsv_file=join(root, 'in', 'subjects.tsv'),
        base_dir=join(working_dir, 'T1VolumeParcellation'),
        parameters=parameters
    )
    pipeline.build()
    pipeline.run(plugin='MultiProc', plugin_args={'n_procs': 4}, bypass_check=True)

    out_files = [join(root, 'out/caps/subjects/sub-ADNI018S4696/ses-M00/t1/spm/dartel/group-UnitTest/atlas_statistics',
                      'sub-ADNI018S4696_ses-M00_T1w_segm-graymatter_space-Ixi549Space_modulated-on_probability_space-'
                      + atlas + '_map-graymatter_statistics.tsv')
                 for atlas in pipeline.parameters['atlases']]
    ref_files = [join(root, 'ref/sub-ADNI018S4696_ses-M00_T1w_segm-graymatter_space-Ixi549Space_modulated-on_probability_space-'
                      + atlas + '_map-graymatter_statistics.tsv')
                 for atlas in pipeline.parameters['atlases']]

    for i in range(len(out_files)):
        out_csv = pds.read_csv(out_files[i], sep='\t')
        ref_csv = pds.read_csv(ref_files[i], sep='\t')
        assert np.allclose(np.array(out_csv.mean_scalar), np.array(ref_csv.mean_scalar), rtol=1e-8, equal_nan=True)

    clean_folder(join(root, 'out', 'caps'), recreate=False)
    clean_folder(join(working_dir, 'T1VolumeParcellation'), recreate=False)


def test_run_DWIPreprocessingUsingT1(cmdopt):
    from os.path import dirname, join, abspath
    from clinica.pipelines.dwi_preprocessing_using_t1.dwi_preprocessing_using_t1_pipeline import DwiPreprocessingUsingT1

    working_dir = cmdopt
    root = dirname(abspath(join(abspath(__file__), pardir)))
    root = join(root, 'data', 'DWIPreprocessingUsingT1')

    # Remove old instance of UT
    clean_folder(join(root, 'out', 'caps'), recreate=True)
    clean_folder(join(working_dir, 'DWIPreprocessingUsingT1'))

    parameters = {'low_bval': 5}
    pipeline = DwiPreprocessingUsingT1(
        bids_directory=join(root, 'in', 'bids'),
        caps_directory=join(root, 'out', 'caps'),
        tsv_file=join(root, 'in', 'subjects.tsv'),
        base_dir=join(working_dir, 'DWIPreprocessingUsingT1'),
        parameters=parameters,
    )
    pipeline.build()
    pipeline.run(plugin='MultiProc', plugin_args={'n_procs': 4}, bypass_check=True)

    # Assert
    out_file = join(root, 'out', 'caps', 'subjects', 'sub-CAPP01001TMM', 'ses-M00', 'dwi', 'preprocessing', 'sub-CAPP01001TMM_ses-M00_dwi_space-T1w_preproc.nii.gz')
    ref_file = join(root, 'ref', 'sub-CAPP01001TMM_ses-M00_dwi_space-T1w_preproc.nii.gz')

    assert similarity_measure(out_file, ref_file, 0.97)

    # Delete out/caps folder
    clean_folder(join(root, 'out', 'caps'), recreate=False)
    clean_folder(join(working_dir, 'DWIPreprocessingUsingT1'), recreate=False)


def test_run_DWIPreprocessingUsingPhaseDiffFieldmap(cmdopt):
    from os.path import dirname, join, abspath
    import warnings
    from clinica.pipelines.dwi_preprocessing_using_phasediff_fieldmap.dwi_preprocessing_using_phasediff_fieldmap_pipeline import DwiPreprocessingUsingPhaseDiffFieldmap
    warnings.filterwarnings("ignore")

    working_dir = cmdopt
    root = dirname(abspath(join(abspath(__file__), pardir)))
    root = join(root, 'data', 'DWIPreprocessingUsingPhaseDiffFieldmap')

    # Remove old instance of UT
    clean_folder(join(root, 'out', 'caps'))
    clean_folder(join(working_dir, 'DWIPreprocessingUsingPhaseDiffFieldmap'))

    parameters = {'low_bval': 5}
    pipeline = DwiPreprocessingUsingPhaseDiffFieldmap(
        bids_directory=join(root, 'in', 'bids'),
        caps_directory=join(root, 'out', 'caps'),
        tsv_file=join(root, 'in', 'subjects.tsv'),
        base_dir=join(working_dir, 'DWIPreprocessingUsingPhaseDiffFieldmap'),
        parameters=parameters)
    pipeline.build()
    pipeline.run(plugin='MultiProc', plugin_args={'n_procs': 4}, bypass_check=True)

    # Assert
    out_file = join(root, 'out', 'caps', 'subjects', 'sub-CAPP01001TMM', 'ses-M00', 'dwi', 'preprocessing', 'sub-CAPP01001TMM_ses-M00_dwi_space-b0_preproc.nii.gz')
    ref_file = join(root, 'ref', 'sub-CAPP01001TMM_ses-M00_dwi_space-b0_preproc.nii.gz')

    assert similarity_measure(out_file, ref_file, 0.95)

    # Delete out/caps folder
    clean_folder(join(root, 'out', 'caps'), recreate=False)
    clean_folder(join(working_dir, 'DWIPreprocessingUsingPhaseDiffFieldmap'), recreate=False)


def test_run_DWIDTI(cmdopt):
    from os.path import dirname, join, abspath
    import shutil
    import pandas as pds
    import numpy as np
    from clinica.pipelines.dwi_dti.dwi_dti_pipeline import DwiDti

    working_dir = cmdopt
    root = dirname(abspath(join(abspath(__file__), pardir)))
    root = join(root, 'data', 'DWIDTI')

    clean_folder(join(root, 'out', 'caps'), recreate=False)
    clean_folder(join(working_dir, 'DWIDTI'))
    shutil.copytree(join(root, 'in', 'caps'), join(root, 'out', 'caps'))

    pipeline = DwiDti(
        caps_directory=join(root, 'out', 'caps'),
        tsv_file=join(root, 'in', 'subjects.tsv'),
        base_dir=join(working_dir, 'DWIDTI')
    )
    pipeline.build()
    pipeline.run(plugin='MultiProc', plugin_args={'n_procs': 4}, bypass_check=True)

    # Check files
    subject_id = 'sub-CAPP01001TMM'
    maps = ['AD', 'FA', 'MD', 'RD']
    out_files = [join(root, 'out', 'caps', 'subjects', subject_id, 'ses-M00', 'dwi', 'dti_based_processing', 'atlas_statistics', subject_id + '_ses-M00_dwi_space-JHUDTI81_res-1x1x1_map-' + m + '_statistics.tsv')
                 for m in maps]
    ref_files = [join(root, 'ref', subject_id + '_ses-M00_dwi_space-JHUDTI81_res-1x1x1_map-' + m + '_statistics.tsv')
                 for m in maps]

    for i in range(len(out_files)):
        out_csv = pds.read_csv(out_files[i], sep='\t')
        out_mean_scalar = np.array(out_csv.mean_scalar)
        ref_csv = pds.read_csv(ref_files[i], sep='\t')
        ref_mean_scalar = np.array(ref_csv.mean_scalar)

        assert np.allclose(out_mean_scalar, ref_mean_scalar, rtol=0.025, equal_nan=True)

    clean_folder(join(root, 'out', 'caps'), recreate=False)
    clean_folder(join(working_dir, 'DWIDTI'), recreate=False)


def test_run_DWIConnectome(cmdopt):
    from os.path import dirname, join, abspath
    import shutil
    from clinica.pipelines.dwi_connectome.dwi_connectome_pipeline import DwiConnectome

    # Initialization
    working_dir = join(abspath(cmdopt), 'DWIConnectome')
    root = dirname(abspath(join(abspath(__file__), pardir)))
    root = join(root, 'data', 'DWIConnectome')
    in_tsv = join(root, 'in', 'subjects.tsv')
    out_caps_dir = join(root, 'out', 'caps')

    subject_id = 'sub-HMTC20110506MEMEPPAT27'
    session_id = 'ses-M00'

    clean_folder(out_caps_dir, recreate=False)
    clean_folder(working_dir)
    shutil.copytree(join(root, 'in', 'caps'), out_caps_dir)

    parameters = {'n_tracks': 1000}
    pipeline = DwiConnectome(
        caps_directory=out_caps_dir,
        tsv_file=in_tsv,
        base_dir=working_dir,
        parameters=parameters
    )
    pipeline.build()
    pipeline.run(plugin='MultiProc', plugin_args={'n_procs': 4}, bypass_check=True)

    # Check files
    atlases = ['desikan', 'destrieux']

    out_fod_file = join(root, 'out', 'caps', 'subjects', subject_id, session_id, 'dwi', 'connectome_based_processing',
                        subject_id + '_' + session_id + '_dwi_space-b0_model-CSD_diffmodel.nii.gz')
    ref_fod_file = join(root, 'ref',
                        subject_id + '_' + session_id + '_dwi_space-b0_model-CSD_diffmodel.nii.gz')

    out_parc_files = [join(root, 'out', 'caps', 'subjects', subject_id, session_id, 'dwi', 'connectome_based_processing',
                           subject_id + '_' + session_id + '_dwi_space-b0_atlas-' + a + '_parcellation.nii.gz')
                      for a in atlases]
    ref_parc_files = [join(root, 'ref',
                           subject_id + '_' + session_id + '_dwi_space-b0_atlas-' + a + '_parcellation.nii.gz')
                      for a in atlases]

    assert similarity_measure(out_fod_file, ref_fod_file, 0.97)

    for i in range(len(out_parc_files)):
        assert similarity_measure(out_parc_files[i], ref_parc_files[i], 0.955)

    clean_folder(join(root, 'out', 'caps'), recreate=False)
    clean_folder(working_dir, recreate=False)


def test_run_fMRIPreprocessing(cmdopt):
    from os.path import dirname, join, abspath
    import shutil
    from clinica.pipelines.fmri_preprocessing.fmri_preprocessing_pipeline import fMRIPreprocessing

    working_dir = cmdopt
    root = dirname(abspath(join(abspath(__file__), pardir)))
    root = join(root, 'data', 'fMRIPreprocessing')

    clean_folder(join(root, 'out', 'caps'), recreate=False)
    clean_folder(join(working_dir, 'fMRIPreprocessing'))
    shutil.copytree(join(root, 'in', 'caps'), join(root, 'out', 'caps'))

    pipeline = fMRIPreprocessing(
        bids_directory=join(root, 'in', 'bids'),
        caps_directory=join(root, 'out', 'caps'),
        tsv_file=join(root, 'in', 'subjects.tsv'),
        base_dir=join(working_dir, 'fMRIPreprocessing'),
    )
    pipeline.build()
    pipeline.run(bypass_check=True)

    subject_id = 'sub-01001TMM'
    out_files = [join(root, 'out', 'caps', 'subjects', subject_id, 'ses-M00', 'fmri', 'preprocessing', subject_id + '_ses-M00_task-rest_bold_space-Ixi549Space_preproc.nii.gz'),
                 join(root, 'out', 'caps', 'subjects', subject_id, 'ses-M00', 'fmri', 'preprocessing', subject_id + '_ses-M00_task-rest_bold_space-meanBOLD_preproc.nii.gz')]

    ref_files = [join(root, 'ref', subject_id + '_ses-M00_task-rest_bold_space-Ixi549Space_preproc.nii.gz'),
                 join(root, 'ref', subject_id + '_ses-M00_task-rest_bold_space-meanBOLD_preproc.nii.gz')]

    for i in range(len(out_files)):
        assert similarity_measure(out_files[i], ref_files[i], 0.99)

    clean_folder(join(root, 'out', 'caps'), recreate=False)
    clean_folder(join(working_dir, 'fMRIPreprocessing'), recreate=False)


def test_run_PETVolume(cmdopt):
    from os.path import dirname, join, abspath
    import shutil
    from clinica.pipelines.pet_volume.pet_volume_pipeline import PETVolume

    working_dir = cmdopt
    root = dirname(abspath(join(abspath(__file__), pardir)))
    root = join(root, 'data', 'PETVolume')

    clean_folder(join(root, 'out', 'caps'), recreate=False)
    clean_folder(join(working_dir, 'PETVolume'))
    shutil.copytree(join(root, 'in', 'caps'), join(root, 'out', 'caps'))

    parameters = {
        'group_id': 'UnitTest'
    }
    pipeline = PETVolume(
        bids_directory=join(root, 'in', 'bids'),
        caps_directory=join(root, 'out', 'caps'),
        tsv_file=join(root, 'in', 'subjects.tsv'),
        base_dir=join(working_dir, 'PETVolume'),
        parameters=parameters
    )
    pipeline.build()
    pipeline.run(plugin='MultiProc', plugin_args={'n_procs': 4}, bypass_check=True)

    subjects = ['sub-ADNI011S4105', 'sub-ADNI023S4020', 'sub-ADNI035S4082', 'sub-ADNI128S4832']
    out_files = [join(root, 'out/caps/subjects/' + sub + '/ses-M00/pet/preprocessing/group-UnitTest',
                      sub + '_ses-M00_task-rest_acq-fdg_pet_space-Ixi549Space_suvr-pons_mask-brain_fwhm-8mm_pet.nii.gz')
                 for sub in subjects]
    ref_files = [join(root, 'ref', sub + '_ses-M00_task-rest_acq-fdg_pet_space-Ixi549Space_suvr-pons_mask-brain_fwhm-8mm_pet.nii.gz')
                 for sub in subjects]

    for i in range(len(out_files)):
        assert likeliness_measure(out_files[i], ref_files[i], (1e-2, 0.25), (1e-1, 0.001))

    clean_folder(join(root, 'out', 'caps'), recreate=False)
    clean_folder(join(working_dir, 'PETVolume'), recreate=False)


def test_run_StatisticsSurface(cmdopt):
    from clinica.pipelines.statistics_surface.statistics_surface_pipeline import StatisticsSurface
    from os.path import dirname, join, abspath
    import shutil
    import numpy as np
    from scipy.io import loadmat

    working_dir = cmdopt
    root = dirname(abspath(join(abspath(__file__), pardir)))
    root = join(root, 'data', 'StatisticsSurface')

    clean_folder(join(root, 'out', 'caps'), recreate=False)
    clean_folder(join(working_dir, 'StatisticsSurface'))
    shutil.copytree(join(root, 'in', 'caps'), join(root, 'out', 'caps'))

    parameters = {
        'design_matrix': '1 + group + age + sex',
        'contrast': 'group',
        'str_format': '%s %s %s %f %s',
        'group_label': 'UnitTest',
        'glm_type': 'group_comparison',
        'custom_file': '@subject/@session/t1/freesurfer_cross_sectional/@subject_@session/surf/@hemi.thickness.fwhm@fwhm.fsaverage.mgh',
        'feature_label': 'cortical_thickness',
        'full_width_at_half_maximum': 20,
        'threshold_uncorrected_pvalue': 0.001,
        'threshold_corrected_pvalue': 0.05,
        'cluster_threshold': 0.001
    }
    pipeline = StatisticsSurface(
        caps_directory=join(root, 'out', 'caps'),
        tsv_file=join(root, 'in', 'subjects.tsv'),
        base_dir=join(working_dir, 'StatisticsSurface'),
        parameters=parameters
    )
    pipeline.build()
    pipeline.run(plugin='MultiProc', plugin_args={'n_procs': 8}, bypass_check=True)

    # Check files
    out_file = join(root, 'out/caps/groups/group-UnitTest/statistics/surfstat_group_comparison/group-UnitTest_AD-lt-CN_measure-cortical_thickness_fwhm-20_correctedPValue.mat')
    ref_file = join(root, 'ref/group-UnitTest_AD-lt-CN_measure-cortical_thickness_fwhm-20_correctedPValue.mat')

    out_file_mat = loadmat(out_file)['correctedpvaluesstruct']
    ref_file_mat = loadmat(ref_file)['correctedpvaluesstruct']
    for i in range(4):
        assert np.allclose(out_file_mat[0][0][i], ref_file_mat[0][0][i], rtol=1e-8, equal_nan=True)
    clean_folder(join(root, 'out', 'caps'), recreate=False)
    clean_folder(join(working_dir, 'StatisticsSurface'), recreate=False)


def test_run_PETSurfaceCrossSectional(cmdopt):
    from os.path import dirname, join, abspath
    import shutil
    import nibabel as nib
    import numpy as np
    from clinica.pipelines.pet_surface.pet_surface_pipeline import PetSurface

    working_dir = cmdopt
    root = dirname(abspath(join(abspath(__file__), pardir)))
    root = join(root, 'data', 'PETSurface')

    clean_folder(join(root, 'out', 'caps'), recreate=False)
    clean_folder(join(working_dir, 'PETSurface'))
    shutil.copytree(join(root, 'in', 'caps'), join(root, 'out', 'caps'))

    parameters = {
        'longitudinal': False
    }
    pipeline = PetSurface(
        bids_directory=join(root, 'in', 'bids'),
        caps_directory=join(root, 'out', 'caps'),
        tsv_file=join(root, 'in', 'subjects.tsv'),
        base_dir=join(working_dir, 'PETSurface'),
        parameters=parameters
    )
    pipeline.build()
    pipeline.run(bypass_check=True)

    # Check files
    out_files = [join(root, 'out/caps/subjects/sub-ADNI011S4105/ses-M00/pet/surface',
                      'sub-ADNI011S4105_ses-M00_task-rest_acq-fdg_pet_space-fsaverage_suvr-pons_pvc-iy_hemi-'
                      + h + '_fwhm-' + str(f) + '_projection.mgh')
                 for h in ['lh', 'rh']
                 for f in [0, 5, 10, 15, 20, 25]]
    ref_files = [join(root, 'ref/sub-ADNI011S4105_ses-M00_task-rest_acq-fdg_pet_space-fsaverage_suvr-pons_pvc-iy_hemi-'
                      + h + '_fwhm-' + str(f) + '_projection.mgh')
                 for h in ['lh', 'rh']
                 for f in [0, 5, 10, 15, 20, 25]]

    for i in range(len(out_files)):
        assert np.allclose(np.squeeze(nib.load(out_files[i]).get_data()),
                           np.squeeze(nib.load(ref_files[i]).get_data()),
                           rtol=3e-2, equal_nan=True)

    clean_folder(join(root, 'out', 'caps'), recreate=False)
    clean_folder(join(working_dir, 'PETSurface'), recreate=False)


# def test_run_PETSurfaceLongitudinal(cmdopt):
#     from os.path import dirname, join, abspath
#     import shutil
#     import nibabel as nib
#     import numpy as np
#     from clinica.pipelines.pet_surface.pet_surface_pipeline import PetSurface
#
#     working_dir = cmdopt
#     root = dirname(abspath(join(abspath(__file__), pardir)))
#     root = join(root, 'data', 'PETSurfaceLongitudinal')
#
#     clean_folder(join(root, 'out', 'caps'), recreate=False)
#     clean_folder(join(working_dir, 'PETSurfaceLongitudinal'))
#     shutil.copytree(join(root, 'in', 'caps'), join(root, 'out', 'caps'))
#
#     parameters = {
#         'longitudinal': True
#     }
#     pipeline = PetSurface(
#         bids_directory=join(root, 'in', 'bids'),
#         caps_directory=join(root, 'out', 'caps'),
#         tsv_file=join(root, 'in', 'subjects.tsv'),
#         base_dir=join(working_dir, 'PETSurfaceLongitudinal'),
#         parameters=parameters
#     )
#     pipeline.build()
#     pipeline.run(bypass_check=True)
#
#     # Check files
#     part_id = 'sub-ADNI041S1260'
#     sess_id = 'ses-M24'
#     long_id = 'long-M00M06M12M18M24'
#     image_id = part_id + '_' + sess_id + '_' + long_id
#     out_files = [join(root, 'out', 'caps', 'subjects', part_id, sess_id, 'pet', long_id, 'surface_longitudinal',
#                       image_id + '_task-rest_acq-fdg_pet_space-fsaverage_suvr-pons_pvc-iy_hemi-'
#                       + h + '_fwhm-' + str(f) + '_projection.mgh')
#                  for h in ['lh', 'rh']
#                  for f in [0, 5, 10, 15, 20, 25]]
#     ref_files = [join(root, 'ref',
#                       image_id + '_task-rest_acq-fdg_pet_space-fsaverage_suvr-pons_pvc-iy_hemi-'
#                       + h + '_fwhm-' + str(f) + '_projection.mgh')
#                  for h in ['lh', 'rh']
#                  for f in [0, 5, 10, 15, 20, 25]]
#
#     # Tolerance values were taken from PETSurface - Cross-sectional case
#     for i in range(len(out_files)):
#         assert np.allclose(np.squeeze(nib.load(out_files[i]).get_data()),
#                            np.squeeze(nib.load(ref_files[i]).get_data()),
#                            rtol=3e-2, equal_nan=True)
#     clean_folder(join(root, 'out', 'caps'), recreate=False)


def test_run_WorkflowsML(cmdopt):
    from clinica.pipelines.machine_learning.ml_workflows import (RegionBasedRepHoldOutLogisticRegression,
                                                                 VertexBasedRepHoldOutDualSVM,
                                                                 RegionBasedRepHoldOutRandomForest,
                                                                 VoxelBasedKFoldDualSVM)
    from os.path import dirname, join, abspath
    import shutil
    import warnings
    warnings.filterwarnings("ignore", category=DeprecationWarning)
    warnings.filterwarnings("ignore", category=UserWarning)
    warnings.filterwarnings("ignore", category=FutureWarning)

    root = dirname(abspath(join(abspath(__file__), pardir)))
    root = join(root, 'data', 'WorkflowsML')
    root_input = dirname(abspath(join(abspath(__file__), pardir)))
    root_input = join(root_input, 'data', 'InputsML')

    caps_dir = join(root_input, 'in', 'caps')
    tsv = join(root_input, 'in', 'subjects.tsv')
    diagnoses_tsv = join(root_input, 'in', 'diagnosis.tsv')
    group_id = 'allADNIdartel'

    output_dir1 = join(root, 'out', 'VertexBasedRepHoldOutDualSVM')
    clean_folder(output_dir1, recreate=True)
    wf1 = VertexBasedRepHoldOutDualSVM(caps_dir, tsv, diagnoses_tsv, group_id, output_dir1, image_type='fdg', fwhm=20,
                                       n_threads=8, n_iterations=10, grid_search_folds=3, test_size=0.3)
    wf1.run()
    shutil.rmtree(output_dir1)

    output_dir2 = join(root, 'out', 'RegionBasedRepHoldOutLogisticRegression')
    clean_folder(output_dir2, recreate=True)
    wf2 = RegionBasedRepHoldOutLogisticRegression(caps_dir, tsv, diagnoses_tsv, group_id, 'fdg', 'AICHA', output_dir2,
                                                  n_threads=8, n_iterations=10, grid_search_folds=3, test_size=0.3)
    wf2.run()
    shutil.rmtree(output_dir2)

    output_dir3 = join(root, 'out', 'RegionBasedRepHoldOutRandomForest')
    clean_folder(output_dir3, recreate=True)
    wf3 = RegionBasedRepHoldOutRandomForest(caps_dir, tsv, diagnoses_tsv, group_id, 'T1', 'AAL2', output_dir3,
                                            n_threads=8, n_iterations=10, grid_search_folds=3, test_size=0.3)
    wf3.run()
    shutil.rmtree(output_dir3)

    output_dir4 = join(root, 'out', 'VoxelBasedKFoldDualSVM')
    clean_folder(output_dir4, recreate=True)
    wf4 = VoxelBasedKFoldDualSVM(caps_dir, tsv, diagnoses_tsv, group_id, 'fdg', output_dir4, fwhm=8, n_threads=8,
                                 n_folds=5, grid_search_folds=3)
    wf4.run()
    shutil.rmtree(output_dir4)


def test_run_SpatialSVM(cmdopt):
    from os.path import dirname, join, abspath
    import shutil
    import numpy as np
    import nibabel as nib
    from clinica.pipelines.machine_learning_spatial_svm.spatial_svm_pipeline import SpatialSVM

    working_dir = cmdopt
    root = dirname(abspath(join(abspath(__file__), pardir)))
    root = join(root, 'data', 'SpatialSVM')

    # Remove potential residual of previous UT
    clean_folder(join(root, 'out', 'caps'), recreate=False)
    clean_folder(join(working_dir, 'SpatialSVM'), recreate=False)

    # Copy necessary data from in to out
    shutil.copytree(join(root, 'in', 'caps'), join(root, 'out', 'caps'))

    parameters = {
        'group_label': 'ADNIbl',
        'orig_input_data': 't1-volume'
    }
    # Instantiate pipeline and run()
    pipeline = SpatialSVM(
        caps_directory=join(root, 'out', 'caps'),
        tsv_file=join(root, 'in', 'subjects.tsv'),
        base_dir=join(working_dir, 'SpatialSVM'),
        parameters=parameters
    )
    pipeline.build()
    pipeline.run(plugin='MultiProc', plugin_args={'n_procs': 4}, bypass_check=True)

    # Check output vs ref
    subjects = ['sub-ADNI011S0023', 'sub-ADNI013S0325']
    out_data_REG_NIFTI = [nib.load(join(root,
                                        'out', 'caps', 'subjects', sub, 'ses-M00',
                                        'machine_learning', 'input_spatial_svm', 'group-ADNIbl',
                                        sub + '_ses-M00_T1w_segm-graymatter_space-Ixi549Space_modulated-on_spatialregularization.nii.gz')).get_data()
                          for sub in subjects]
    ref_data_REG_NIFTI = [nib.load(join(root, 'ref', sub + '_ses-M00_T1w_segm-graymatter_space-Ixi549Space_modulated-on_spatialregularization.nii.gz')).get_data()
                          for sub in subjects]
    for i in range(len(out_data_REG_NIFTI)):
        assert np.allclose(out_data_REG_NIFTI[i], ref_data_REG_NIFTI[i],
                           rtol=1e-3, equal_nan=True)

    # Remove data in out folder
    clean_folder(join(root, 'out', 'caps'), recreate=True)
    clean_folder(join(working_dir, 'SpatialSVM'), recreate=False)


def test_run_T1Linear(cmdopt):
    from os.path import dirname, join, abspath
    import shutil
    from clinica.pipelines.t1_linear.t1_linear_pipeline import T1Linear
    import nibabel as nib
    import numpy as np

    working_dir = cmdopt
    root = dirname(abspath(join(abspath(__file__), pardir)))
    root = join(root, 'data', 'T1Linear')

    # Remove potential residual of previous UT
    clean_folder(join(working_dir, 'T1Linear'))
    clean_folder(join(root, 'out', 'caps'), recreate=False)
    shutil.copytree(join(root, 'in', 'caps'), join(root, 'out', 'caps'))

    parameters = {
        'crop_image': True
    }
    # Instantiate pipeline
    pipeline = T1Linear(
        bids_directory=join(root, 'in', 'bids'),
        caps_directory=join(root, 'out', 'caps'),
        tsv_file=join(root, 'in', 'subjects.tsv'),
        base_dir=join(working_dir, 'T1Linear'),
        parameters=parameters
    )
    pipeline.run(plugin='MultiProc', plugin_args={'n_procs': 4}, bypass_check=True)

    # Check output vs ref

    out_folder = join(root, 'out')
    ref_folder = join(root, 'out')

    compare_folders(out_folder, ref_folder, shared_folder_name='caps')

    clean_folder(join(root, 'out', 'caps'), recreate=False)
    clean_folder(join(working_dir, 'T1Linear'), recreate=False)


def test_run_DLPrepareData(cmdopt):
    from os.path import dirname, join, abspath
    import shutil
    from clinica.pipelines.deeplearning_prepare_data.deeplearning_prepare_data_pipeline import DeepLearningPrepareData
    import nibabel as nib
    import numpy as np

    working_dir = cmdopt
    root = dirname(abspath(join(abspath(__file__), pardir)))
    root = join(root, 'data', 'DeepLearningPrepareData')

    # Remove potential residual of previous UT
    clean_folder(join(working_dir, 'DeepLearningPrepareData'))
    clean_folder(join(root, 'out', 'caps'), recreate=False)

    # Copy necessary data from in to out
    shutil.copytree(join(root, 'in', 'caps'), join(root, 'out', 'caps'))

    # Test the transformation of the complete T1 MRI
    parameters = {
        'extract_method': 'image'
    }
    # Instantiate pipeline
    pipeline = DeepLearningPrepareData(
        caps_directory=join(root, 'out', 'caps'),
        tsv_file=join(root, 'in', 'subjects.tsv'),
        base_dir=join(working_dir, 'DeepLearningPrepareData'),
        parameters=parameters
    )
    pipeline.run(plugin='MultiProc', plugin_args={'n_procs': 4}, bypass_check=True)

    # Test the patch extraction
    parameters = {
        'extract_method': 'patch',
        'patch_size': 50,
        'stride_size': 50
    }
    # Instantiate pipeline
    pipeline = DeepLearningPrepareData(
        caps_directory=join(root, 'out', 'caps'),
        tsv_file=join(root, 'in', 'subjects.tsv'),
        base_dir=join(working_dir, 'DeepLearningPrepareData'),
        parameters=parameters
    )
    pipeline.run(plugin='MultiProc', plugin_args={'n_procs': 4}, bypass_check=True)

    # Test the slice extraction
    parameters = {
        'extract_method': 'slice',
        'slice_mode': 'rgb',
        'slice_direction': 0
    }
    # Instantiate pipeline
    pipeline = DeepLearningPrepareData(
        caps_directory=join(root, 'out', 'caps'),
        tsv_file=join(root, 'in', 'subjects.tsv'),
        base_dir=join(working_dir, 'DeepLearningPrepareData'),
        parameters=parameters
    )
    pipeline.run(plugin='MultiProc', plugin_args={'n_procs': 4}, bypass_check=True)
    # Check output vs ref

    out_folder = join(root, 'out')
    ref_folder = join(root, 'out')

    compare_folders(out_folder, ref_folder, shared_folder_name='caps')

    clean_folder(join(root, 'out', 'caps'), recreate=False)
    clean_folder(join(working_dir, 'DeepLearningPrepareData'), recreate=False)



def test_run_StatisticsVolume(cmdopt):
    from os.path import dirname, join, abspath
    import shutil
    import numpy as np
    import nibabel as nib
    from clinica.pipelines.statistics_volume.statistics_volume_pipeline import StatisticsVolume

    working_dir = cmdopt
    root = dirname(abspath(join(abspath(__file__), pardir)))
    root = join(root, 'data', 'StatisticsVolume')

    # Remove potential residual of previous UT
    clean_folder(join(root, 'out', 'caps'), recreate=False)
    clean_folder(join(working_dir, 'StatisticsVolume'), recreate=False)

    # Copy necessary data from in to out
    shutil.copytree(join(root, 'in', 'caps'), join(root, 'out', 'caps'))

    # Instantiate pipeline and run()
    parameters = {
        'contrast': 'group',
        'feature_type': 'fdg',
        'group_id': 'UnitTest',
        'cluster_threshold': 0.001,
        'group_id_caps': None,
        'full_width_at_half_maximum': 8
    }

    pipeline = StatisticsVolume(
        caps_directory=join(root, 'out', 'caps'),
        tsv_file=join(root, 'in', 'group-UnitTest_covariates.tsv'),
        base_dir=join(working_dir, 'StatisticsVolume'),
        parameters=parameters
    )

    pipeline.run(plugin='MultiProc', plugin_args={'n_procs': 2}, bypass_check=True)

    output_t_stat = join(root, 'out',
                         'caps', 'groups', 'group-UnitTest', 'statistics_volume', 'group_comparison_measure-fdg',
                         'group-UnitTest_CN-lt-AD_measure-fdg_fwhm-8_TStatistics.nii')
    ref_t_stat = join(root, 'ref',
                      'caps', 'groups', 'group-UnitTest', 'statistics_volume', 'group_comparison_measure-fdg',
                      'group-UnitTest_CN-lt-AD_measure-fdg_fwhm-8_TStatistics.nii')

    assert np.allclose(nib.load(output_t_stat).get_data(),
                       nib.load(ref_t_stat).get_data())

    # Remove data in out folder
    clean_folder(join(root, 'out', 'caps'), recreate=True)
    clean_folder(join(working_dir, 'StatisticsVolume'), recreate=False)


def test_run_StatisticsVolumeCorrection(cmdopt):
    from clinica.pipelines.statistics_volume_correction.statistics_volume_correction_pipeline import StatisticsVolumeCorrection
    from os.path import dirname, join, abspath
    import shutil

    working_dir = cmdopt
    root = dirname(abspath(join(abspath(__file__), pardir)))
    root = join(root, 'data', 'StatisticsVolumeCorrection')

    # Remove potential residual of previous UT
    clean_folder(join(root, 'out', 'caps'), recreate=False)
    clean_folder(join(working_dir, 'StatisticsVolumeCorrection'), recreate=False)

    # Copy necessary data from in to out
    shutil.copytree(join(root, 'in', 'caps'), join(root, 'out', 'caps'))

    # Instantiate pipeline and run()
    parameters = {
        't_map': 'group-UnitTest_AD-lt-CN_measure-fdg_fwhm-8_TStatistics.nii',
        'height_threshold': 3.2422,
        'FWEp': 4.928,
        'FDRp': 4.693,
        'FWEc': 206987,
        'FDRc': 206987,
        'n_cuts': 15
    }
    pipeline = StatisticsVolumeCorrection(
        caps_directory=join(root, 'out', 'caps'),
        base_dir=join(working_dir, 'StatisticsVolumeCorrection'),
        parameters=parameters
    )
    pipeline.build()
    pipeline.run(plugin='MultiProc', plugin_args={'n_procs': 4}, bypass_check=True)
    compare_folders(join(root, 'out'), join(root, 'ref'), 'caps')

    # Remove data in out folder
    clean_folder(join(root, 'out', 'caps'), recreate=True)
    clean_folder(join(working_dir, 'StatisticsVolumeCorrection'), recreate=False)


# def test_run_T1FreeSurferTemplate(cmdopt):
#     # Data for this functional test comes from https://openneuro.org/datasets/ds000204
#     # sub-01 was duplicated into to sub-02 with one session in order to test the "one time point" case
#     import shutil
#     from os.path import dirname, join, abspath
#     from clinica.pipelines.t1_freesurfer_longitudinal.t1_freesurfer_template_pipeline import T1FreeSurferTemplate
#
#     working_dir = cmdopt
#     root = dirname(abspath(join(abspath(__file__), pardir)))
#     root = join(root, 'data', 'T1FreeSurferTemplate')
#
#     # Remove potential residual of previous tests
#     clean_folder(join(root, 'out', 'caps'), recreate=False)
#     clean_folder(join(working_dir, 'T1FreeSurferTemplate'))
#
#     # Copy necessary data from in to out
#     shutil.copytree(join(root, 'in', 'caps'), join(root, 'out', 'caps'))
#
#     pipeline = T1FreeSurferTemplate(
#         caps_directory=join(root, 'out', 'caps'),
#         tsv_file=join(root, 'in', 'subjects.tsv'),
#         base_dir=join(working_dir, 'T1FreeSurferTemplate'),
#     )
#     pipeline.base_dir = join(working_dir, 'T1FreeSurferTemplate')
#     pipeline.run(plugin='MultiProc',
#                  plugin_args={'n_procs': 2},
#                  bypass_check=True)
#
#     # We only check that folders are the same meaning that FreeSurfer finished without error
#     # surf/ folder is ignored because it contains sym links that makes hard to check with ref data
#     # (sym links of ref data are ignored after rsync on CI machines)
#     def path_to_caps_fs(part_id, long_id):
#         import os
#         output_folder = os.path.join('caps', 'subjects', part_id, long_id, 'freesurfer_unbiased_template')
#         return output_folder
#
#     for (p_id, l_id) in zip(['sub-01', 'sub-02'], ['long-20112015', 'long-2011']):
#         compare_folders(join(root, 'out'), join(root, 'ref'),
#                         join(path_to_caps_fs(p_id, l_id), p_id + '_' + l_id, 'label'))
#         compare_folders(join(root, 'out'), join(root, 'ref'),
#                         join(path_to_caps_fs(p_id, l_id), p_id + '_' + l_id, 'mri'))
#         compare_folders(join(root, 'out'), join(root, 'ref'),
#                         join(path_to_caps_fs(p_id, l_id), p_id + '_' + l_id, 'stats'))
#
#     clean_folder(join(root, 'out', 'caps'), recreate=False)
#     clean_folder(join(working_dir, 'T1FreeSurferTemplate'), recreate=False)


# def test_run_T1FreeSurferLongitudinalCorrection(cmdopt):
#     # Data for this functional test comes from https://openneuro.org/datasets/ds000204
#     import shutil
#     from os.path import dirname, join, abspath
#     from clinica.pipelines.t1_freesurfer_longitudinal.t1_freesurfer_longitudinal_correction_pipeline import T1FreeSurferLongitudinalCorrection
#
#     working_dir = cmdopt
#     root = dirname(abspath(join(abspath(__file__), pardir)))
#     root = join(root, 'data', 'T1FreeSurferLongitudinalCorrection')
#
#     # Remove potential residual of previous tests
#     clean_folder(join(root, 'out', 'caps'), recreate=False)
#     clean_folder(join(working_dir, 'T1FreeSurferLongitudinalCorrection'))
#
#     # Copy necessary data from in to out
#     shutil.copytree(join(root, 'in', 'caps'), join(root, 'out', 'caps'))
#
#     pipeline = T1FreeSurferLongitudinalCorrection(
#         caps_directory=join(root, 'out', 'caps'),
#         tsv_file=join(root, 'in', 'subjects.tsv'),
#         base_dir=join(working_dir, 'T1FreeSurferLongitudinalCorrection'),
#     )
#     pipeline.base_dir = join(working_dir, 'T1FreeSurferLongitudinalCorrection')
#     pipeline.run(bypass_check=True)
#
#     # We only check that folders are the same meaning that FreeSurfer finished without error
#     # surf/ folder is ignored because it contains sym links that makes hard to check with ref data
#     # (sym links of ref data are ignored after rsync on CI machines)
#     def path_to_caps_fs(part_id, sess_id, long_id):
#         import os
#         output_folder = os.path.join('caps', 'subjects', part_id, sess_id, 't1', long_id, 'freesurfer_longitudinal')
#         return output_folder
#
#     compare_folders(join(root, 'out'), join(root, 'ref'),
#                     join(path_to_caps_fs('sub-01', 'ses-2011', 'long-20112015'),
#                          'regional_measures'))
#     compare_folders(join(root, 'out'), join(root, 'ref'),
#                     join(path_to_caps_fs('sub-01', 'ses-2011', 'long-20112015'),
#                          'sub-01_ses-2011.long.sub-01_long-20112015', 'label'))
#     compare_folders(join(root, 'out'), join(root, 'ref'),
#                     join(path_to_caps_fs('sub-01', 'ses-2011', 'long-20112015'),
#                          'sub-01_ses-2011.long.sub-01_long-20112015', 'mri'))
#     compare_folders(join(root, 'out'), join(root, 'ref'),
#                     join(path_to_caps_fs('sub-01', 'ses-2011', 'long-20112015'),
#                          'sub-01_ses-2011.long.sub-01_long-20112015', 'stats'))
#
#     clean_folder(join(root, 'out', 'caps'), recreate=False)
#     clean_folder(join(working_dir, 'T1FreeSurferLongitudinalCorrection'), recreate=False)
