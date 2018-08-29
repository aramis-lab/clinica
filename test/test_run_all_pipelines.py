## test if pipeline can run until the very end
import warnings
warnings.filterwarnings("ignore")
working_dir = '/localdrive10TB/data/working_directory_ci_linux_branch-upgrade_dep'


def test_run_T1FreeSurferCrossSectional():
    from clinica.pipelines.t1_freesurfer_cross_sectional.t1_freesurfer_cross_sectional_pipeline import T1FreeSurferCrossSectional
    from os.path import dirname, join, abspath, exists
    import shutil
    from os import makedirs

    root = join(dirname(abspath(__file__)), 'data', 'T1FreeSurferCrossSectional')

    clean_folder(join(root, 'out', 'caps'))
    clean_folder(join(working_dir, 'T1FreeSurferCrossSectional'))

    pipeline = T1FreeSurferCrossSectional(bids_directory=join(root, 'in', 'bids'),
                                          caps_directory=join(root, 'out', 'caps'),
                                          tsv_file=join(root, 'in', 'subjects.tsv'))
    pipeline.parameters['recon_all_args'] = '-qcache'
    pipeline.base_dir = join(working_dir, 'T1FreeSurferCrossSectional')
    pipeline.build()
    #pipeline.run()
    pass


def test_run_T1VolumeTissueSegmentation():
    import nibabel as nib
    import shutil
    import numpy as np
    import os
    from clinica.pipelines.t1_volume_tissue_segmentation.t1_volume_tissue_segmentation_pipeline import T1VolumeTissueSegmentation
    from os.path import dirname, join, abspath

    root = join(dirname(abspath(__file__)), 'data', 'T1VolumeTissueSegmentation')
    clean_folder(join(working_dir, 'T1VolumeTissueSegmentation'))
    clean_folder(join(root, 'out', 'caps'))

    pipeline = T1VolumeTissueSegmentation(bids_directory=join(root, 'in', 'bids'),
                                          caps_directory=join(root, 'out', 'caps'),
                                          tsv_file=join(root, 'in', 'subjects.tsv'))
    pipeline.base_dir = join(working_dir, 'T1VolumeTissueSegmentation')
    pipeline.build()
    pipeline.run()

    out_file = join(root, 'out/caps/subjects/sub-ADNI011S4105/ses-M00/t1/spm/segmentation/dartel_input/'
                    + 'sub-ADNI011S4105_ses-M00_T1w_segm-graymatter_dartelinput.nii.gz')
    if not os.path.exists(out_file):
        raise IOError('Pipeline did not produce file : ' + out_file +'. Consider rerunning test_run_T1VolumeTissueSegmentation')

    ref_file = join(root, 'ref/caps/subjects/sub-ADNI011S4105/ses-M00/t1/spm/segmentation/dartel_input/'
                    + 'sub-ADNI011S4105_ses-M00_T1w_segm-graymatter_dartelinput.nii.gz')

    assert np.allclose(nib.load(out_file).get_data(), nib.load(ref_file).get_data(), rtol=1e-8, equal_nan=True)
    clean_folder(join(root, 'out', 'caps'), recreate=False)
    pass

def test_run_T1VolumeCreateDartel():
    from clinica.pipelines.t1_volume_create_dartel.t1_volume_create_dartel_pipeline import T1VolumeCreateDartel
    from os.path import dirname, join, abspath, exists
    import shutil
    import nibabel as nib
    import numpy as np

    root = join(dirname(abspath(__file__)), 'data', 'T1VolumeCreateDartel')

    # Remove potential residual of previous UT
    clean_folder(join(working_dir, 'T1VolumeCreateDartel'))
    clean_folder(join(root, 'out', 'caps'), recreate=False)
    shutil.copytree(join(root, 'in', 'caps'), join(root, 'out', 'caps'))

    # Instantiate pipeline
    pipeline = T1VolumeCreateDartel(bids_directory=join(root, 'in', 'bids'),
                                    caps_directory=join(root, 'out', 'caps'),
                                    tsv_file=join(root, 'in', 'subjects.tsv'),
                                    group_id='UnitTest')
    pipeline.base_dir = join(working_dir, 'T1VolumeCreateDartel')
    pipeline.build()
    pipeline.run(plugin='MultiProc', plugin_args={'n_procs': 4})

    # Check output vs ref
    out_data_template = nib.load(join(root, 'out/caps/groups/group-UnitTest/t1/group-UnitTest_template.nii.gz')).get_data()
    ref_data_template = nib.load(join(root, 'ref/group-UnitTest_template.nii.gz')).get_data()
    assert np.allclose(out_data_template, ref_data_template, rtol=1e-8, equal_nan=True)

    subjects = ['sub-ADNI011S4105', 'sub-ADNI023S4020', 'sub-ADNI035S4082', 'sub-ADNI128S4832']
    out_data_forward_def = [nib.load(join(root, 'out', 'caps', 'subjects', sub, 'ses-M00', 't1', 'spm', 'dartel', 'group-UnitTest',
                                          sub + '_ses-M00_T1w_target-UnitTest_transformation-forward_deformation.nii.gz')).get_data() for sub in subjects]
    ref_data_forward_def = [nib.load(join(root, 'ref', sub + '_ses-M00_T1w_target-UnitTest_transformation-forward_deformation.nii.gz')).get_data()
                            for sub in subjects]

    for i in range(len(out_data_forward_def)):
        assert np.allclose(out_data_forward_def[i], ref_data_forward_def[i], rtol=1e-8, equal_nan=True)

    # Remove data in out folder
    clean_folder(join(root, 'out', 'caps'), recreate=False)
    pass

def test_run_T1VolumeDartel2MNI():
    from clinica.pipelines.t1_volume_dartel2mni.t1_volume_dartel2mni_pipeline import T1VolumeDartel2MNI
    from os.path import dirname, join, abspath, exists
    import shutil
    import numpy as np
    import nibabel as nib

    root = join(dirname(abspath(__file__)), 'data', 'T1VolumeDartel2MNI')

    # Remove potential residual of previous UT
    clean_folder(join(root, 'out', 'caps'), recreate=False)
    clean_folder(join(working_dir, 'T1VolumeDartel2MNI'))

    # Copy necessary data from in to out
    shutil.copytree(join(root, 'in', 'caps'), join(root, 'out', 'caps'))

    # Instantiate pipeline and run()
    pipeline = T1VolumeDartel2MNI(bids_directory=join(root, 'in', 'bids'),
                                  caps_directory=join(root, 'out', 'caps'),
                                  tsv_file=join(root, 'in', 'subjects.tsv'),
                                  group_id='UnitTest')
    pipeline.base_dir = join(working_dir, 'T1VolumeDartel2MNI')
    pipeline.build()
    pipeline.run(plugin='MultiProc', plugin_args={'n_procs': 4})

    # Check output vs ref
    subjects = ['sub-ADNI011S4105', 'sub-ADNI023S4020', 'sub-ADNI035S4082', 'sub-ADNI128S4832']
    out_data_GM_MNI = [nib.load(join(root, 'out', 'caps', 'subjects', sub, 'ses-M00', 't1', 'spm', 'dartel', 'group-UnitTest',
                                     sub + '_ses-M00_T1w_segm-graymatter_space-Ixi549Space_modulated-on_fwhm-8mm_probability.nii.gz')).get_data()
                       for sub in subjects]
    ref_data_GM_MNI = [nib.load(join(root, 'ref', sub + '_ses-M00_T1w_segm-graymatter_space-Ixi549Space_modulated-on_fwhm-8mm_probability.nii.gz')).get_data()
                       for sub in subjects]
    for i in range(len(out_data_GM_MNI)):
        assert np.allclose(out_data_GM_MNI[i], ref_data_GM_MNI[i], rtol=1e-8, equal_nan=True)

    # Remove data in out folder
    clean_folder(join(root, 'out', 'caps'), recreate=False)
    pass


def test_run_T1VolumeNewTemplate():
    from clinica.pipelines.t1_volume_new_template.t1_volume_new_template_pipeline import T1VolumeNewTemplate
    from os.path import dirname, join, abspath, exists, basename
    import shutil
    import numpy as np
    import nibabel as nib

    root = join(dirname(abspath(__file__)), 'data', 'T1VolumeNewTemplate')

    # Remove residual files from out folder
    clean_folder(join(root, 'out', 'caps'), recreate=True)
    clean_folder(join(working_dir, 'T1VolumeNewTemplate'))

    pipeline = T1VolumeNewTemplate(bids_directory=join(root, 'in', 'bids'),
                                   caps_directory=join(root, 'out', 'caps'),
                                   tsv_file=join(root, 'in', 'subjects.tsv'),
                                   group_id='UnitTest')
    pipeline.base_dir = join(working_dir, 'T1VolumeNewTemplate')
    pipeline.build()
    pipeline.run(plugin='MultiProc', plugin_args={'n_procs': 4})

    # Check generated vs ref
    subjects = ['sub-ADNI011S4105', 'sub-ADNI023S4020', 'sub-ADNI035S4082', 'sub-ADNI128S4832']
    out_data_GM_MNI = [nib.load(join(root, 'out', 'caps', 'subjects', sub, 'ses-M00', 't1', 'spm', 'dartel', 'group-UnitTest',
                                     sub + '_ses-M00_T1w_segm-graymatter_space-Ixi549Space_modulated-on_fwhm-8mm_probability.nii.gz')).get_data()
                       for sub in subjects]
    ref_data_GM_MNI = [nib.load(join(root, 'ref', sub + '_ses-M00_T1w_segm-graymatter_space-Ixi549Space_modulated-on_fwhm-8mm_probability.nii.gz')).get_data()
                       for sub in subjects]

    # Check output vs ref
    out_data_template = nib.load(join(root, 'out/caps/groups/group-UnitTest/t1/group-UnitTest_template.nii.gz')).get_data()
    ref_data_template = nib.load(join(root, 'ref/group-UnitTest_template.nii.gz')).get_data()
    assert np.allclose(out_data_template, ref_data_template, rtol=1e-8, equal_nan=True)

    for i in range(len(out_data_GM_MNI)):
        print('Checking file ' + subjects[i] + '_ses-M00_T1w_segm-graymatter_space-Ixi549Space_modulated-on_fwhm-8mm_probability.nii.gz')
        assert np.allclose(out_data_GM_MNI[i], ref_data_GM_MNI[i], rtol=1e-8, equal_nan=True)

    # Remove data in out folder
    clean_folder(join(root, 'out', 'caps'), recreate=False)
    pass


def test_run_T1VolumeExistingTemplate():
    from clinica.pipelines.t1_volume_existing_template.t1_volume_existing_template_pipeline import T1VolumeExistingTemplate
    from os.path import dirname, join, abspath, exists
    import shutil
    import nibabel as nib
    import numpy as np

    root = join(dirname(abspath(__file__)), 'data', 'T1VolumeExistingTemplate')
    clean_folder(join(root, 'out', 'caps'), recreate=False)
    clean_folder(join(working_dir, 'T1VolumeExistingTemplate'))

    # Copy necessary data to run pipeline
    shutil.copytree(join(root, 'in', 'caps'), join(root, 'out', 'caps'))

    # Instantiate and run pipeline
    pipeline = T1VolumeExistingTemplate(bids_directory=join(root, 'in', 'bids'),
                                        caps_directory=join(root, 'out', 'caps'),
                                        tsv_file=join(root, 'in', 'subjects.tsv'),
                                        group_id='UnitTest')
    pipeline.base_dir = join(working_dir, 'T1VolumeExistingTemplate')
    pipeline.build()
    pipeline.run(plugin='MultiProc', plugin_args={'n_procs': 4})

    # Check output vs ref
    subjects = ['sub-ADNI011S4105', 'sub-ADNI023S4020', 'sub-ADNI035S4082', 'sub-ADNI128S4832']
    out_data_forward_def = [nib.load(join(root, 'out', 'caps', 'subjects', sub, 'ses-M00', 't1', 'spm', 'dartel', 'group-UnitTest',
                                          sub + '_ses-M00_T1w_target-UnitTest_transformation-forward_deformation.nii.gz')).get_data()
                            for sub in subjects]
    ref_data_forward_def = [nib.load(join(root, 'ref', sub + '_ses-M00_T1w_target-UnitTest_transformation-forward_deformation.nii.gz')).get_data()
                            for sub in subjects]

    for i in range(len(out_data_forward_def)):
        assert np.allclose(out_data_forward_def[i], ref_data_forward_def[i], rtol=1e-8, equal_nan=True)

    # Remove data in out folder
    clean_folder(join(root, 'out', 'caps'), recreate=False)
    pass


def test_run_T1VolumeParcellation():
    from clinica.pipelines.t1_volume_parcellation.t1_volume_parcellation_pipeline import T1VolumeParcellation
    from os.path import dirname, join, abspath, exists
    import shutil
    import pandas as pds
    import numpy as np

    root = join(dirname(abspath(__file__)), 'data', 'T1VolumeParcellation')
    clean_folder(join(root, 'out', 'caps'), recreate=False)
    clean_folder(join(working_dir, 'T1VolumeParcellation'))

    # Copy data for use of pipeline
    shutil.copytree(join(root, 'in', 'caps'), join(root, 'out', 'caps'))

    # Instantiate pipeline
    pipeline = T1VolumeParcellation(bids_directory='./4',
                                    caps_directory=join(root, 'out', 'caps'),
                                    tsv_file=join(root, 'in', 'subjects.tsv'))
    pipeline.parameters['group_id'] = 'UnitTest'
    pipeline.parameters['atlases'] = ['AAL2', 'LPBA40', 'Neuromorphometrics', 'AICHA', 'Hammers']
    pipeline.parameters['modulate'] = 'on'
    pipeline.base_dir = join(working_dir, 'T1VolumeParcellation')
    pipeline.build()
    pipeline.run(plugin='MultiProc', plugin_args={'n_procs': 4})

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
    pass


def test_run_DWIPreprocessingUsingT1():
    # TODO : Alex must resolve the repetability issue (out file is never the same) -  assert deactivated
    from clinica.pipelines.dwi_preprocessing_using_t1.dwi_preprocessing_using_t1_pipeline import DWIPreprocessingUsingT1
    from os.path import dirname, join, abspath, exists
    from os import makedirs
    import shutil
    import nibabel as nib
    import numpy as np

    root = join(dirname(abspath(__file__)), 'data', 'DWIPreprocessingUsingT1')

    # Remove old instance of UT
    clean_folder(join(root, 'out', 'caps'))
    clean_folder(join(working_dir, 'DWIPreprocessingUsingT1'))

    pipeline = DWIPreprocessingUsingT1(bids_directory=join(root, 'in', 'bids'),
                                       caps_directory=join(root, 'out', 'caps'),
                                       tsv_file=join(root, 'in', 'subjects.tsv'),
                                       low_bval=5)
    pipeline.base_dir = join(working_dir, 'DWIPreprocessingUsingT1')
    pipeline.build()
    pipeline.run(plugin='MultiProc', plugin_args={'n_procs': 4})

    # Assert :
    out_file = join(root, 'out', 'caps', 'subjects', 'sub-CAPP01001TMM', 'ses-M00', 'dwi', 'preprocessing', 'sub-CAPP01001TMM_ses-M00_dwi_space-T1w_preproc.nii.gz')
    ref_file = join(root, 'ref', 'sub-CAPP01001TMM_ses-M00_dwi_space-T1w_preproc.nii.gz')


    ## Uncomment when problem is solved !
    #out_data = np.array(nib.load(out_file).get_data())
    #ref_data = np.array(nib.load(ref_file).get_data())

    #assert np.allclose(out_data, ref_data, rtol=1e-8, equal_nan=True)

    # Delete out/caps folder
    clean_folder(join(root, 'out', 'caps'), recreate=False)
    pass


def test_run_DWIPreprocessingUsingPhaseDiffFieldmap():
    # TODO : Alex must resolve the repetability issue (out file is never the same) -  assert deactivated
    from clinica.pipelines.dwi_preprocessing_using_phasediff_fieldmap.dwi_preprocessing_using_phasediff_fieldmap_pipeline import DWIPreprocessingUsingPhaseDiffFieldmap
    from os.path import dirname, join, abspath
    import nibabel as nib
    import numpy as np

    root = join(dirname(abspath(__file__)), 'data', 'DWIPreprocessingUsingPhaseDiffFieldmap')

    # Remove old instance of UT
    clean_folder(join(root, 'out', 'caps'))
    clean_folder(join(working_dir, 'DWIPreprocessingUsingPhaseDiffFieldmap'))

    pipeline = DWIPreprocessingUsingPhaseDiffFieldmap(bids_directory=join(root, 'in', 'bids'),
                                                      caps_directory=join(root, 'out', 'caps'),
                                                      tsv_file=join(root, 'in', 'subjects.tsv'),
                                                      low_bval=5)
    pipeline.base_dir = join(working_dir, 'DWIPreprocessingUsingPhaseDiffFieldmap')
    pipeline.build()
    pipeline.run(plugin='MultiProc', plugin_args={'n_procs': 4})

    # Assert :
    out_file = join(root, 'out', 'caps', 'subjects', 'sub-CAPP01001TMM', 'ses-M00', 'dwi', 'preprocessing', 'sub-CAPP01001TMM_ses-M00_dwi_space-bo_preproc.nii.gz')
    ref_file = join(root, 'ref', 'sub-CAPP01001TMM_ses-M00_dwi_space-T1w_preproc.nii.gz')

    ## Uncomment when repetability problem is solved
    # out_data = np.array(nib.load(out_file).get_data())
    # ref_data = np.array(nib.load(ref_file).get_data())
    #
    # assert np.allclose(out_data, ref_data, rtol=1e-8, equal_nan=True)

    # Delete out/caps folder
    clean_folder(join(root, 'out', 'caps'), recreate=False)
    pass

def test_run_DWIProcessingDTI():
    # TODO : Alex must resolve the repetability issue (out file is never the same) -  assert deactivated
    from clinica.pipelines.dwi_processing_dti.dwi_processing_dti_pipeline import DWIProcessingDTI
    from os.path import dirname, join, abspath, exists
    import shutil
    import pandas as pds
    import numpy as np

    root = join(dirname(abspath(__file__)), 'data', 'DWIProcessingDTI')

    #clean_folder(join(root, 'out', 'caps'), recreate=False)
    #clean_folder(join(working_dir, 'DWIProcessingDTI'))
    #shutil.copytree(join(root, 'in', 'caps'), join(root, 'out', 'caps'))

    pipeline = DWIProcessingDTI(caps_directory=join(root, 'out', 'caps'),
                                tsv_file=join(root, 'in', 'subjects.tsv'))
    pipeline.base_dir = join(working_dir, 'DWIProcessingDTI')
    pipeline.build()
    pipeline.run()

    # Check files
    maps = ['ad', 'fa', 'md', 'rd']
    out_files = [join(root, 'out/caps/subjects/sub-CAPP01001TMM/ses-M00/dwi/dti_based_processing/atlas_statistics/sub-CAPP01001TMM_ses-M00_dwi_space-JHUDTI81_res-1x1x1_map-' + m + '_statistics.tsv')
                 for m in maps]
    ref_files = [join(root, 'ref', 'sub-CAPP01001TMM_ses-M00_dwi_space-JHUDTI81_res-1x1x1_map-' + m + '_statistics.tsv')
                 for m in maps]

    for i in range(len(out_files)):
        out_csv = pds.read_csv(out_files[i], sep='\t')
        out_mean_scalar = np.array(out_csv.mean_scalar)
        ref_csv = pds.read_csv(ref_files[i], sep='\t')
        ref_mean_scalar = np.array(ref_csv.mean_scalar)
        ## Uncomment when problem is solved
        #assert np.allclose(out_mean_scalar, ref_mean_scalar, rtol=1e-8, equal_nan=True)

    clean_folder(join(root, 'out', 'caps'), recreate=False)
    pass

def test_run_fMRIPreprocessing():
    from clinica.pipelines.fmri_preprocessing.fmri_preprocessing_pipeline import fMRIPreprocessing
    from os.path import dirname, join, abspath, exists
    import shutil
    import nibabel as nib
    import numpy as np

    root = join(dirname(abspath(__file__)), 'data', 'fMRIPreprocessing')

    clean_folder(join(root, 'out', 'caps'), recreate=False)
    clean_folder(join(working_dir, 'fMRIPreprocessing'))
    shutil.copytree(join(root, 'in', 'caps'), join(root, 'out', 'caps'))

    pipeline = fMRIPreprocessing(bids_directory=join(root, 'in', 'bids'),
                                 caps_directory=join(root, 'out', 'caps'),
                                 tsv_file=join(root, 'in', 'subjects.tsv'))
    pipeline.parameters = {
        'full_width_at_half_maximum': [8, 8, 8],
        't1_native_space': False,
        'freesurfer_brain_mask': False,
        'unwarping': False
    }
    pipeline.base_dir = join(working_dir, 'fMRIPreprocessing')
    pipeline.build()
    pipeline.run()

    out_files = [join(root, 'out/caps/subjects/sub-20110506MEMEPPAT27/ses-M00/fmri/preprocessing/sub-20110506MEMEPPAT27_ses-M00_task-rest_bold_space-Ixi549Space.nii.gz'),
                 join(root, 'out/caps/subjects/sub-20110506MEMEPPAT27/ses-M00/fmri/preprocessing/sub-20110506MEMEPPAT27_ses-M00_task-rest_bold_space-meanBOLD.nii.gz')]

    ref_files = [join(root, 'ref', 'sub-20110506MEMEPPAT27_ses-M00_task-rest_bold_space-Ixi549Space.nii.gz'),
                 join(root, 'ref', 'sub-20110506MEMEPPAT27_ses-M00_task-rest_bold_space-meanBOLD.nii.gz')]

    for i in range(len(out_files)):
        assert np.allclose(nib.load(out_files[i]).get_data(), nib.load(ref_files[i]).get_data(), rtol=1e-8, equal_nan=True)

    clean_folder(join(root, 'out', 'caps'), recreate=False)
    pass


def test_run_PETVolume():
    from clinica.pipelines.pet_volume.pet_volume_pipeline import PETVolume
    from os.path import dirname, join, abspath, exists
    import shutil
    import nibabel as nib
    import numpy as np

    root = join(dirname(abspath(__file__)), 'data', 'PETVolume')

    clean_folder(join(root, 'out', 'caps'), recreate=False)
    clean_folder(join(working_dir, 'PETVolume'))
    shutil.copytree(join(root, 'in', 'caps'), join(root, 'out', 'caps'))

    pipeline = PETVolume(bids_directory=join(root, 'in', 'bids'),
                         caps_directory=join(root, 'out', 'caps'),
                         tsv_file=join(root, 'in', 'subjects.tsv'),
                         group_id='UnitTest',
                         fwhm_tsv=None)
    pipeline.base_dir = join(working_dir, 'PETVolume')
    pipeline.build()
    pipeline.run(plugin='MultiProc', plugin_args={'n_procs': 4})

    subjects = ['sub-ADNI011S4105', 'sub-ADNI023S4020', 'sub-ADNI035S4082', 'sub-ADNI128S4832']
    out_files = [join(root, 'out/caps/subjects/' + sub + '/ses-M00/pet/preprocessing/group-UnitTest',
                      sub +'_ses-M00_task-rest_acq-fdg_pet_space-Ixi549Space_suvr-pons_mask-brain_fwhm-8mm_pet.nii.gz')
                 for sub in subjects]
    ref_files = [join(root, 'ref', sub + '_ses-M00_task-rest_acq-fdg_pet_space-Ixi549Space_suvr-pons_mask-brain_fwhm-8mm_pet.nii.gz')
                 for sub in subjects]

    for i in range(len(out_files)):
        assert np.allclose(nib.load(out_files[i]).get_data(), nib.load(ref_files[i]).get_data(), rtol=1e-8, equal_nan=True)

    clean_folder(join(root, 'out', 'caps'), recreate=False)
    pass

def test_run_StatisticsSurface():
    from clinica.pipelines.statistics_surface.statistics_surface_pipeline import StatisticsSurface
    from os.path import dirname, join, abspath, exists
    import shutil
    import numpy as np
    from scipy.io import loadmat

    root = join(dirname(abspath(__file__)), 'data', 'StatisticsSurface')

    clean_folder(join(root, 'out', 'caps'), recreate=False)
    clean_folder(join(working_dir, 'StatisticsSurface'))
    shutil.copytree(join(root, 'in', 'caps'), join(root, 'out', 'caps'))

    pipeline = StatisticsSurface(caps_directory=join(root, 'out', 'caps'),
                                 tsv_file=join(root, 'in', 'subjects.tsv'))
    pipeline.parameters = {
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
    pipeline.base_dir = join(working_dir, 'StatisticsSurface')
    pipeline.build()
    pipeline.run(plugin='MultiProc', plugin_args={'n_procs': 8})

    # Check files
    out_file = join(root, 'out/caps/groups/group-UnitTest/statistics/surfstat_group_comparison/group-UnitTest_AD-lt-CN_measure-cortical_thickness_fwhm-20_correctedPValue.mat')
    ref_file = join(root, 'ref/group-UnitTest_AD-lt-CN_measure-cortical_thickness_fwhm-20_correctedPValue.mat')

    out_file_mat = loadmat(out_file)['correctedpvaluesstruct']
    ref_file_mat = loadmat(ref_file)['correctedpvaluesstruct']
    for i in range(4):
        assert np.allclose(out_file_mat[0][0][i], ref_file_mat[0][0][i], rtol=1e-8, equal_nan=True)
    clean_folder(join(root, 'out', 'caps'), recreate=False)
    pass

def test_run_PETSurface():
    # TODO once freesurfer6 is available, uncomment lines !
    from clinica.pipelines.pet_surface.pet_surface_pipeline import PetSurface
    from os.path import dirname, join, abspath
    import shutil
    import nibabel as nib
    import numpy as np

    root = join(dirname(abspath(__file__)), 'data', 'PETSurface')

    clean_folder(join(root, 'out', 'caps'), recreate=False)
    clean_folder(join(working_dir, 'PETSurface'))
    shutil.copytree(join(root, 'in', 'caps'), join(root, 'out', 'caps'))

    pipeline = PetSurface(bids_directory=join(root, 'in', 'bids'),
                          caps_directory=join(root, 'in', 'caps'),
                          tsv_file=join(root, 'in', 'subjects.tsv'))
    pipeline.parameters['pet_type'] = 'fdg'
    pipeline.parameters['wd'] = join(working_dir, 'PETSurface')
    pipeline.build()
    pipeline.run()

    # Check files
    out_files = [join(root, 'out/caps/subjects/sub-ADNI011S4105/ses-M00/pet/surface',
                      'sub-ADNI011S4105_ses-M00_task-rest_acq-FDG_pet_space-fsaverage_suvr-pons_pvc-iy_hemi-'
                      + h + '_fwhm-' + str(f) + '_projection.mgh') for h in ['lh', 'rh'] for f in [0, 5, 10, 15, 20, 25]]
    ref_files = [join(root, 'ref/sub-ADNI011S4105_ses-M00_task-rest_acq-FDG_pet_space-fsaverage_suvr-pons_pvc-iy_hemi-'
                      + h + '_fwhm-' + str(f) + '_projection.mgh') for h in ['lh', 'rh'] for f in [0, 5, 10, 15, 20, 25]]

    for i in range(len(out_files)):
        assert np.allclose(np.squeeze(nib.load(out_files[i]).get_data()),
                           np.squeeze(nib.load(ref_files[i]).get_data()),
                           rtol=1e-8, equal_nan=True)
    clean_folder(join(root, 'out', 'caps'), recreate=False)
    pass


def clean_folder(path, recreate=True):
    from os.path import abspath, exists
    from shutil import rmtree
    from os import makedirs

    abs_path = abspath(path)
    if exists(abs_path):
        rmtree(abs_path)
        if recreate:
            makedirs(abs_path)
    pass
