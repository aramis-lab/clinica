# coding: utf8

# Small unit tests for all pipelines
##
# test if instantiation and building of workflows is working
import warnings
warnings.filterwarnings("ignore")


def test_instantiate_T1FreeSurferCrossSectional():
    from clinica.pipelines.t1_freesurfer_cross_sectional.t1_freesurfer_cross_sectional_pipeline import T1FreeSurferCrossSectional
    from os.path import dirname, join, abspath

    root = dirname(abspath(__file__))
    pipeline = T1FreeSurferCrossSectional(bids_directory=join(root, 'data', 'T1FreeSurferCrossSectional', 'in', 'bids'),
                                          caps_directory=join(root, 'data', 'T1FreeSurferCrossSectional', 'in', 'caps'),
                                          tsv_file=join(root, 'data', 'T1FreeSurferCrossSectional', 'in', 'subjects.tsv'))
    pipeline.parameters['recon_all_args'] = '-qcache'
    pipeline.build()



def test_instantiate_T1VolumeTissueSegmentation():
    from clinica.pipelines.t1_volume_tissue_segmentation.t1_volume_tissue_segmentation_pipeline import T1VolumeTissueSegmentation
    from os.path import dirname, join, abspath

    root = dirname(abspath(__file__))
    pipeline = T1VolumeTissueSegmentation(bids_directory=join(root, 'data', 'T1VolumeTissueSegmentation', 'in', 'bids'),
                                          caps_directory=join(root, 'data', 'T1VolumeTissueSegmentation', 'in', 'caps'),
                                          tsv_file=join(root, 'data', 'T1VolumeTissueSegmentation', 'in', 'subjects.tsv'))
    pipeline.build()



def test_instantiate_T1VolumeCreateDartel():
    from clinica.pipelines.t1_volume_create_dartel.t1_volume_create_dartel_pipeline import T1VolumeCreateDartel
    from os.path import dirname, join, abspath

    root = dirname(abspath(__file__))
    pipeline = T1VolumeCreateDartel(bids_directory=join(root, 'data', 'T1VolumeCreateDartel', 'in', 'bids'),
                                    caps_directory=join(root, 'data', 'T1VolumeCreateDartel', 'in', 'caps'),
                                    tsv_file=join(root, 'data', 'T1VolumeCreateDartel', 'in', 'subjects.tsv'),
                                    group_id='UnitTest')
    pipeline.build()


def test_instantiate_T1VolumeDartel2MNI():
    from clinica.pipelines.t1_volume_dartel2mni.t1_volume_dartel2mni_pipeline import T1VolumeDartel2MNI
    from os.path import dirname, join, abspath

    root = dirname(abspath(__file__))
    pipeline = T1VolumeDartel2MNI(bids_directory=join(root, 'data', 'T1VolumeDartel2MNI', 'in', 'bids'),
                                  caps_directory=join(root, 'data', 'T1VolumeDartel2MNI', 'in', 'caps'),
                                  tsv_file=join(root, 'data', 'T1VolumeDartel2MNI', 'in', 'subjects.tsv'),
                                  group_id='UnitTest')
    pipeline.build()



def test_instantiate_T1VolumeNewTemplate():
    from clinica.pipelines.t1_volume_new_template.t1_volume_new_template_pipeline import T1VolumeNewTemplate
    from os.path import dirname, join, abspath

    root = dirname(abspath(__file__))
    pipeline = T1VolumeNewTemplate(bids_directory=join(root, 'data', 'T1VolumeNewTemplate', 'in', 'bids'),
                                   caps_directory=join(root, 'data', 'T1VolumeNewTemplate', 'in', 'caps'),
                                   tsv_file=join(root, 'data', 'T1VolumeNewTemplate', 'in', 'subjects.tsv'),
                                   group_id='UnitTest')
    pipeline.build()


def test_instantiate_T1VolumeExistingDartel():
    from clinica.pipelines.t1_volume_existing_dartel.t1_volume_existing_dartel_pipeline import T1VolumeExistingDartel
    from os.path import dirname, join, abspath

    root = dirname(abspath(__file__))
    pipeline = T1VolumeExistingDartel(bids_directory=join(root, 'data', 'T1VolumeExistingDartel', 'in', 'bids'),
                                      caps_directory=join(root, 'data', 'T1VolumeExistingDartel', 'in', 'caps'),
                                      tsv_file=join(root, 'data', 'T1VolumeExistingDartel', 'in', 'subjects.tsv'),
                                      group_id='UnitTest')
    pipeline.build()



def test_instantiate_T1VolumeExistingTemplate():
    from clinica.pipelines.t1_volume_existing_template.t1_volume_existing_template_pipeline import T1VolumeExistingTemplate
    from os.path import dirname, join, abspath

    root = join(dirname(abspath(__file__)), 'data', 'T1VolumeExistingTemplate')
    pipeline = T1VolumeExistingTemplate(bids_directory=join(root, 'in', 'bids'),
                                        caps_directory=join(root, 'in', 'caps'),
                                        tsv_file=join(root, 'in', 'subjects.tsv'),
                                        group_id='UnitTest')
    pipeline.build()
    pass


def test_instantiate_T1VolumeParcellation():
    from clinica.pipelines.t1_volume_parcellation.t1_volume_parcellation_pipeline import T1VolumeParcellation
    from os.path import dirname, join, abspath

    root = dirname(abspath(__file__))
    pipeline = T1VolumeParcellation(caps_directory=join(root, 'data', 'T1VolumeParcellation', 'in', 'caps'),
                                    tsv_file=join(root, 'data', 'T1VolumeParcellation', 'in', 'subjects.tsv'))
    pipeline.parameters['group_id'] = 'UnitTest'
    pipeline.parameters['atlases'] = ['AAL2', 'LPBA40', 'Neuromorphometrics', 'AICHA', 'Hammers']
    pipeline.parameters['modulate'] = 'on'
    pipeline.build()



def test_instantiate_DWIPreprocessingUsingT1():
    from clinica.pipelines.dwi_preprocessing_using_t1.dwi_preprocessing_using_t1_pipeline import DwiPreprocessingUsingT1
    from os.path import dirname, join, abspath

    root = dirname(abspath(__file__))
    pipeline = DwiPreprocessingUsingT1(bids_directory=join(root, 'data', 'DWIPreprocessingUsingT1', 'in', 'bids'),
                                       caps_directory=join(root, 'data', 'DWIPreprocessingUsingT1', 'in', 'caps'),
                                       tsv_file=join(root, 'data', 'DWIPreprocessingUsingT1', 'in', 'subjects.tsv'),
                                       low_bval=5)
    pipeline.parameters = {
            'epi_param': dict([('readout_time', 0.14),  ('enc_dir', 'y')]),
    }
    pipeline.build()



def test_instantiate_DWIPreprocessingUsingPhaseDiffFieldmap():
    from clinica.pipelines.dwi_preprocessing_using_phasediff_fieldmap.dwi_preprocessing_using_phasediff_fieldmap_pipeline import DwiPreprocessingUsingPhaseDiffFieldmap
    from os.path import dirname, join, abspath

    root = dirname(abspath(__file__))
    pipeline = DwiPreprocessingUsingPhaseDiffFieldmap(bids_directory=join(root, 'data', 'DWIPreprocessingUsingPhaseDiffFieldmap', 'in', 'bids'),
                                                      caps_directory=join(root, 'data', 'DWIPreprocessingUsingPhaseDiffFieldmap', 'in', 'caps'),
                                                      tsv_file=join(root, 'data', 'DWIPreprocessingUsingPhaseDiffFieldmap', 'in', 'subjects.tsv'),
                                                      low_bval=5)
    pipeline.build()



def test_instantiate_DWIDTI():
    from clinica.pipelines.dwi_dti.dwi_dti_pipeline import DwiDti
    from os.path import dirname, join, abspath

    root = dirname(abspath(__file__))
    pipeline = DwiDti(caps_directory=join(root, 'data', 'DWIDTI', 'in', 'caps'),
                      tsv_file=join(root, 'data', 'DWIDTI', 'in', 'subjects.tsv'))
    pipeline.build()
    pass


def test_instantiate_DWIConnectome():
    from clinica.pipelines.dwi_connectome.dwi_connectome_pipeline import DwiConnectome
    from os.path import dirname, join, abspath

    root = dirname(abspath(__file__))
    pipeline = DwiConnectome(caps_directory=join(root, 'data', 'DWIConnectome', 'in', 'caps'),
                             tsv_file=join(root, 'data', 'DWIConnectome', 'in', 'subjects.tsv'))
    pipeline.parameters = {
        'n_tracks' : 1000
    }
    pipeline.build()



def test_instantiate_fMRIPreprocessing():
    # Need to add json file in BIDS
    from clinica.pipelines.fmri_preprocessing.fmri_preprocessing_pipeline import fMRIPreprocessing
    from os.path import dirname, join, abspath

    root = dirname(abspath(__file__))
    pipeline = fMRIPreprocessing(bids_directory=join(root, 'data', 'fMRIPreprocessing', 'in', 'bids'),
                                 caps_directory=join(root, 'data', 'fMRIPreprocessing', 'in', 'caps'),
                                 tsv_file=join(root, 'data', 'fMRIPreprocessing', 'in', 'subjects.tsv'))
    pipeline.parameters = {
            'full_width_at_half_maximum' : [8, 8, 8],
            't1_native_space'            : True,
            'freesurfer_brain_mask'      : True,
            'unwarping'                  : True
    }
    pipeline.build()



def test_instantiate_PETVolume():
    from clinica.pipelines.pet_volume.pet_volume_pipeline import PETVolume
    from os.path import dirname, join, abspath

    root= dirname(abspath(__file__))
    pipeline = PETVolume(bids_directory=join(root, 'data', 'PETVolume', 'in', 'bids'),
                         caps_directory=join(root, 'data', 'PETVolume', 'in', 'caps'),
                         tsv_file=join(root, 'data', 'PETVolume', 'in', 'subjects.tsv'),
                         group_id='UnitTest',
                         fwhm_tsv=None)
    pipeline.build()



def test_instantiate_StatisticsSurface():
    from clinica.pipelines.statistics_surface.statistics_surface_pipeline import StatisticsSurface
    from os.path import dirname, join, abspath

    root= dirname(abspath(__file__))
    pipeline = StatisticsSurface(caps_directory=join(root, 'data', 'StatisticsSurface', 'in', 'caps'),
                                 tsv_file=join(root, 'data', 'StatisticsSurface', 'in', 'subjects.tsv'))
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
    pipeline.build()



def test_instantiate_PETSurface(tmpdir):
    from clinica.pipelines.pet_surface.pet_surface_pipeline import PetSurface
    from os.path import dirname, join, abspath

    root= dirname(abspath(__file__))
    pipeline = PetSurface(bids_directory=join(root, 'data', 'PETSurface', 'in', 'bids'),
                          caps_directory=join(root, 'data', 'PETSurface', 'in', 'caps'),
                          tsv_file=join(root, 'data', 'PETSurface', 'in', 'subjects.tsv'))
    pipeline.parameters['pet_type'] = 'fdg'
    pipeline.parameters['wd'] = str(tmpdir)
    pipeline.build()


def test_instantiate_InputsML():
    from clinica.pipelines.machine_learning.input import CAPSVoxelBasedInput, CAPSRegionBasedInput, CAPSVertexBasedInput
    from os.path import dirname, join, abspath, exists

    root = join(dirname(abspath(__file__)), 'data', 'InputsML')

    caps_dir = join(root, 'in', 'caps')
    tsv = join(root, 'in', 'subjects.tsv')
    diagnoses_tsv = join(root, 'in', 'diagnosis.tsv')
    group_id = 'allADNIdartel'
    image_type = ['T1', 'fdg']
    atlases = ['AAL2', 'Neuromorphometrics', 'AICHA', 'LPBA40', 'Hammers']
    possible_psf = [0, 5, 10, 15, 20, 25]

    voxel_input = [CAPSVoxelBasedInput(caps_dir, tsv, diagnoses_tsv, group_id, im, fwhm=8)
                   for im in image_type]
    region_input = [CAPSRegionBasedInput(caps_dir, tsv, diagnoses_tsv, group_id, im, at)
                    for im in image_type
                    for at in atlases]
    vertex_input = [CAPSVertexBasedInput(caps_dir, tsv, diagnoses_tsv, group_id, fwhm, 'fdg')
                    for fwhm in possible_psf]

    # Check that each file exists
    for inputs in voxel_input + region_input + vertex_input:
        for file in inputs.get_images():
            if isinstance(file, str):
                assert exists(file)
            elif isinstance(file, list) and len(file) == 2:
                assert exists(file[0])
                assert exists(file[1])
            else:
                raise ValueError('An error occured...')


def test_instantiate_SpatialSVM():
    from clinica.pipelines.machine_learning_spatial_svm.spatial_svm_pipeline import SpatialSVM
    from os.path import dirname, join, abspath

    root = join(dirname(abspath(__file__)), 'data', 'SpatialSVM')
    pipeline = SpatialSVM(caps_directory=join(root, 'in', 'caps'),
                          tsv_file=join(root, 'in', 'subjects.tsv'))
    pipeline.parameters['group_id'] = 'ADNIbl'
    pipeline.parameters['fwhm'] = 4
    pipeline.parameters['image_type'] = 't1'
    pipeline.parameters['pet_type'] = 'fdg'
    pipeline.parameters['no_pvc'] = 'True'
    pipeline.build()
