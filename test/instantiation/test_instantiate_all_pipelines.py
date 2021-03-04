# coding: utf8

import warnings

# Small unit tests for all pipelines
##
# test if instantiation and building of workflows is working
from os import pardir

warnings.filterwarnings("ignore")


def test_instantiate_T1FreeSurferCrossSectional():
    from os.path import abspath, dirname, join

    from clinica.pipelines.t1_freesurfer.t1_freesurfer_pipeline import T1FreeSurfer

    root = dirname(abspath(join(abspath(__file__), pardir)))
    root = join(root, "data", "T1FreeSurfer")

    parameters = {"recon_all_args": "-qcache", "skip_question": False}

    pipeline = T1FreeSurfer(
        bids_directory=join(root, "in", "bids"),
        caps_directory=join(root, "in", "caps"),
        tsv_file=join(root, "in", "subjects.tsv"),
        parameters=parameters,
    )
    pipeline.build()


def test_instantiate_T1VolumeTissueSegmentation():
    from os.path import abspath, dirname, join

    from clinica.pipelines.t1_volume_tissue_segmentation.t1_volume_tissue_segmentation_pipeline import (
        T1VolumeTissueSegmentation,
    )

    root = dirname(abspath(join(abspath(__file__), pardir)))
    root = join(root, "data", "T1VolumeTissueSegmentation")

    pipeline = T1VolumeTissueSegmentation(
        bids_directory=join(root, "in", "bids"),
        caps_directory=join(root, "in", "caps"),
        tsv_file=join(root, "in", "subjects.tsv"),
    )
    pipeline.build()


def test_instantiate_T1VolumeCreateDartel():
    from os.path import abspath, dirname, join

    from clinica.pipelines.t1_volume_create_dartel.t1_volume_create_dartel_pipeline import (
        T1VolumeCreateDartel,
    )

    root = dirname(abspath(join(abspath(__file__), pardir)))
    root = join(root, "data", "T1VolumeCreateDartel")

    parameters = {"group_label": "UnitTest"}
    pipeline = T1VolumeCreateDartel(
        bids_directory=join(root, "in", "bids"),
        caps_directory=join(root, "in", "caps"),
        tsv_file=join(root, "in", "subjects.tsv"),
        parameters=parameters,
    )
    pipeline.build()


def test_instantiate_T1VolumeDartel2MNI():
    from os.path import abspath, dirname, join

    from clinica.pipelines.t1_volume_dartel2mni.t1_volume_dartel2mni_pipeline import (
        T1VolumeDartel2MNI,
    )

    root = dirname(abspath(join(abspath(__file__), pardir)))
    root = join(root, "data", "T1VolumeDartel2MNI")

    parameters = {"group_label": "UnitTest"}
    pipeline = T1VolumeDartel2MNI(
        bids_directory=join(root, "in", "bids"),
        caps_directory=join(root, "in", "caps"),
        tsv_file=join(root, "in", "subjects.tsv"),
        parameters=parameters,
    )
    pipeline.build()


def test_instantiate_T1VolumeRegisterDartel():
    from os.path import abspath, dirname, join

    from clinica.pipelines.t1_volume_register_dartel.t1_volume_register_dartel_pipeline import (
        T1VolumeRegisterDartel,
    )

    root = dirname(abspath(join(abspath(__file__), pardir)))
    root = join(root, "data", "T1VolumeExistingDartel")

    parameters = {"group_label": "UnitTest"}
    pipeline = T1VolumeRegisterDartel(
        bids_directory=join(root, "in", "bids"),
        caps_directory=join(root, "in", "caps"),
        tsv_file=join(root, "in", "subjects.tsv"),
        parameters=parameters,
    )
    pipeline.build()


def test_instantiate_T1VolumeParcellation():
    from os.path import abspath, dirname, join

    from clinica.pipelines.t1_volume_parcellation.t1_volume_parcellation_pipeline import (
        T1VolumeParcellation,
    )

    root = dirname(abspath(join(abspath(__file__), pardir)))
    root = join(root, "data", "T1VolumeParcellation")

    parameters = {"group_label": "UnitTest"}
    pipeline = T1VolumeParcellation(
        caps_directory=join(root, "in", "caps"),
        tsv_file=join(root, "in", "subjects.tsv"),
        parameters=parameters,
    )
    pipeline.build()


def test_instantiate_DWIPreprocessingUsingT1():
    from os.path import abspath, dirname, join

    from clinica.pipelines.dwi_preprocessing_using_t1.dwi_preprocessing_using_t1_pipeline import (
        DwiPreprocessingUsingT1,
    )

    root = dirname(abspath(join(abspath(__file__), pardir)))
    root = join(root, "data", "DWIPreprocessingUsingT1")

    parameters = {"low_bval": 5}
    pipeline = DwiPreprocessingUsingT1(
        bids_directory=join(root, "in", "bids"),
        caps_directory=join(root, "in", "caps"),
        tsv_file=join(root, "in", "subjects.tsv"),
        parameters=parameters,
    )
    pipeline.build()


def test_instantiate_DWIPreprocessingUsingPhaseDiffFieldmap():
    from os.path import abspath, dirname, join

    from clinica.pipelines.dwi_preprocessing_using_fmap.dwi_preprocessing_using_phasediff_fmap_pipeline import (
        DwiPreprocessingUsingPhaseDiffFMap,
    )

    root = dirname(abspath(join(abspath(__file__), pardir)))
    root = join(root, "data", "DWIPreprocessingUsingPhaseDiffFieldmap")

    parameters = {"low_bval": 5}
    pipeline = DwiPreprocessingUsingPhaseDiffFMap(
        bids_directory=join(root, "in", "bids"),
        caps_directory=join(root, "in", "caps"),
        tsv_file=join(root, "in", "subjects.tsv"),
        parameters=parameters,
    )
    pipeline.build()


def test_instantiate_DWIDTI():
    from os.path import abspath, dirname, join

    from clinica.pipelines.dwi_dti.dwi_dti_pipeline import DwiDti

    root = dirname(abspath(join(abspath(__file__), pardir)))
    root = join(root, "data", "DWIDTI")

    pipeline = DwiDti(
        caps_directory=join(root, "in", "caps"),
        tsv_file=join(root, "in", "subjects.tsv"),
    )
    pipeline.build()


def test_instantiate_DWIConnectome():
    from os.path import abspath, dirname, join

    from clinica.pipelines.dwi_connectome.dwi_connectome_pipeline import DwiConnectome

    root = dirname(abspath(join(abspath(__file__), pardir)))
    root = join(root, "data", "DWIConnectome")

    parameters = {"n_tracks": 1000}
    pipeline = DwiConnectome(
        caps_directory=join(root, "in", "caps"),
        tsv_file=join(root, "in", "subjects.tsv"),
        parameters=parameters,
    )
    pipeline.build()


def test_instantiate_PETVolume():
    from os.path import abspath, dirname, join

    from clinica.pipelines.pet_volume.pet_volume_pipeline import PETVolume

    root = dirname(abspath(join(abspath(__file__), pardir)))
    root = join(root, "data", "PETVolume")

    parameters = {
        "group_label": "UnitTest",
        "acq_label": "fdg",
        "suvr_reference_region": "pons",
    }
    pipeline = PETVolume(
        bids_directory=join(root, "in", "bids"),
        caps_directory=join(root, "in", "caps"),
        tsv_file=join(root, "in", "subjects.tsv"),
        parameters=parameters,
    )
    pipeline.build()


def test_instantiate_PETLinear():
    from os.path import abspath, dirname, join

    from clinica.pipelines.pet_linear.pet_linear_pipeline import PETLinear

    root = dirname(abspath(join(abspath(__file__), pardir)))
    root = join(root, "data", "PETLinear")

    parameters = {"acq_label": "fdg", "suvr_reference_region": "cerebellumPons2"}
    pipeline = PETLinear(
        bids_directory=join(root, "in", "bids"),
        caps_directory=join(root, "in", "caps"),
        tsv_file=join(root, "in", "subjects.tsv"),
        parameters=parameters,
    )
    pipeline.build


def test_instantiate_StatisticsSurface():
    from os.path import abspath, dirname, join

    from clinica.pipelines.statistics_surface.statistics_surface_pipeline import (
        StatisticsSurface,
    )

    root = dirname(abspath(join(abspath(__file__), pardir)))
    root = join(root, "data", "StatisticsSurface")

    parameters = {
        # Clinica compulsory parameters
        "group_label": "UnitTest",
        "orig_input_data": "t1-freesurfer",
        "glm_type": "group_comparison",
        "contrast": "group",
        # Optional parameters
        "covariates": "age sex",
    }
    pipeline = StatisticsSurface(
        caps_directory=join(root, "in", "caps"),
        tsv_file=join(root, "in", "subjects.tsv"),
        parameters=parameters,
    )

    pipeline.build()


def test_instantiate_PETSurfaceCrossSectional():
    from os.path import abspath, dirname, join

    from clinica.pipelines.pet_surface.pet_surface_pipeline import PetSurface

    root = dirname(abspath(join(abspath(__file__), pardir)))
    root = join(root, "data", "PETSurface")
    parameters = {
        "acq_label": "FDG",
        "suvr_reference_region": "pons",
        "pvc_psf_tsv": join(root, "in", "subjects.tsv"),
        "longitudinal": False,
    }
    pipeline = PetSurface(
        bids_directory=join(root, "in", "bids"),
        caps_directory=join(root, "in", "caps"),
        tsv_file=join(root, "in", "subjects.tsv"),
        parameters=parameters,
    )
    pipeline.build()


# def test_instantiate_PETSurfaceLongitudinal():
#     from os.path import dirname, join, abspath
#     from clinica.pipelines.pet_surface.pet_surface_pipeline import PetSurface
#
#     root = dirname(abspath(join(abspath(__file__), pardir)))
#     root = join(root, 'data', 'PETSurfaceLongitudinal')
#     parameters = {
#         'acq_label': 'FDG',
#         'suvr_reference_region': 'pons',
#         'pvc_psf_tsv': join(root, 'in', 'subjects.tsv'),
#         'longitudinal': True
#     }
#     pipeline = PetSurface(
#         bids_directory=join(root, 'in', 'bids'),
#         caps_directory=join(root, 'in', 'caps'),
#         tsv_file=join(root, 'in', 'subjects.tsv'),
#         parameters=parameters,
#     )
#     pipeline.build()


def test_instantiate_InputsML():
    from os.path import abspath, dirname, exists, join

    from clinica.pipelines.machine_learning.input import (
        CAPSRegionBasedInput,
        CAPSVertexBasedInput,
        CAPSVoxelBasedInput,
    )

    root = dirname(abspath(join(abspath(__file__), pardir)))
    root = join(root, "data", "InputsML")
    caps_dir = join(root, "in", "caps")
    tsv = join(root, "in", "subjects.tsv")
    diagnoses_tsv = join(root, "in", "diagnosis.tsv")
    group_label = "allADNIdartel"
    image_type = ["T1w", "PET"]
    atlases = ["AAL2", "Neuromorphometrics", "AICHA", "LPBA40", "Hammers"]
    possible_fwhm = [0, 5, 10, 15, 20, 25]

    voxel_input = [
        CAPSVoxelBasedInput(
            {
                "caps_directory": caps_dir,
                "subjects_visits_tsv": tsv,
                "diagnoses_tsv": diagnoses_tsv,
                "group_label": group_label,
                "image_type": im,
                "fwhm": 8,
                "acq_label": "fdg",
                "suvr_reference_region": "pons",
                "use_pvc_data": False,
            }
        )
        for im in image_type
    ]

    region_input = [
        CAPSRegionBasedInput(
            {
                "caps_directory": caps_dir,
                "subjects_visits_tsv": tsv,
                "diagnoses_tsv": diagnoses_tsv,
                "group_label": group_label,
                "image_type": im,
                "atlas": at,
                "acq_label": "fdg",
                "suvr_reference_region": "pons",
                "use_pvc_data": False,
            }
        )
        for im in image_type
        for at in atlases
    ]

    vertex_input = [
        CAPSVertexBasedInput(
            {
                "caps_directory": caps_dir,
                "subjects_visits_tsv": tsv,
                "diagnoses_tsv": diagnoses_tsv,
                "group_label": group_label,
                "image_type": "PET",
                "fwhm": fwhm,
                "acq_label": "fdg",
                "suvr_reference_region": "pons",
            }
        )
        for fwhm in possible_fwhm
    ]

    # Check that each file exists
    for inputs in voxel_input + region_input + vertex_input:
        for file in inputs.get_images():
            if isinstance(file, str):
                assert exists(file)
            elif isinstance(file, list) and len(file) == 2:
                assert exists(file[0])
                assert exists(file[1])
            else:
                raise ValueError("An error occurred...")


def test_instantiate_SpatialSVM():
    from os.path import abspath, dirname, join

    from clinica.pipelines.machine_learning_spatial_svm.spatial_svm_pipeline import (
        SpatialSVM,
    )

    root = dirname(abspath(join(abspath(__file__), pardir)))
    root = join(root, "data", "SpatialSVM")

    parameters = {"group_label": "ADNIbl", "orig_input_data": "t1-volume"}
    pipeline = SpatialSVM(
        caps_directory=join(root, "in", "caps"),
        tsv_file=join(root, "in", "subjects.tsv"),
        parameters=parameters,
    )
    pipeline.build()


def test_instantiate_T1FreeSurferTemplate():
    from os.path import abspath, dirname, join

    from clinica.pipelines.t1_freesurfer_longitudinal.t1_freesurfer_template_pipeline import (
        T1FreeSurferTemplate,
    )

    root = dirname(abspath(join(abspath(__file__), pardir)))
    root = join(root, "data", "T1FreeSurferTemplate")

    pipeline = T1FreeSurferTemplate(
        caps_directory=join(root, "in", "caps"),
        tsv_file=join(root, "in", "subjects.tsv"),
    )
    pipeline.build()


def test_instantiate_T1FreeSurferLongitudinalCorrection():
    from os.path import abspath, dirname, join

    from clinica.pipelines.t1_freesurfer_longitudinal.t1_freesurfer_longitudinal_correction_pipeline import (
        T1FreeSurferLongitudinalCorrection,
    )

    root = dirname(abspath(join(abspath(__file__), pardir)))
    root = join(root, "data", "T1FreeSurferLongitudinalCorrection")

    pipeline = T1FreeSurferLongitudinalCorrection(
        caps_directory=join(root, "in", "caps"),
        tsv_file=join(root, "in", "subjects.tsv"),
    )
    pipeline.build()


def test_instantiate_T1Linear():
    from os.path import abspath, dirname, join

    from clinica.pipelines.t1_linear.t1_linear_pipeline import T1Linear

    root = dirname(abspath(join(abspath(__file__), pardir)))
    root = join(root, "data", "T1Linear")

    parameters = {"uncropped_image": False}

    pipeline = T1Linear(
        bids_directory=join(root, "in", "bids"),
        caps_directory=join(root, "in", "caps"),
        tsv_file=join(root, "in", "subjects.tsv"),
    )
    pipeline.build()


def test_instantiate_DLPrepareData():
    from os.path import abspath, dirname, join

    from clinica.pipelines.deeplearning_prepare_data.deeplearning_prepare_data_pipeline import (
        DeepLearningPrepareData,
    )

    root = dirname(abspath(join(abspath(__file__), pardir)))
    root = join(root, "data", "DeepLearningPrepareData")

    parameters = {"modality": "t1-linear", "extract_method": "image"}
    pipeline = DeepLearningPrepareData(
        caps_directory=join(root, "in", "caps"),
        tsv_file=join(root, "in", "subjects.tsv"),
        parameters=parameters,
    )
    pipeline.build()


def test_instantiate_StatisticsVolume():
    from os.path import abspath, dirname, join

    from clinica.pipelines.statistics_volume.statistics_volume_pipeline import (
        StatisticsVolume,
    )

    root = dirname(abspath(join(abspath(__file__), pardir)))
    root = join(root, "data", "StatisticsVolume")

    # Instantiate pipeline and run()
    parameters = {
        # Clinica compulsory parameters
        "group_label": "UnitTest",
        "orig_input_data": "pet-volume",
        "contrast": "group",
        # Optional arguments for inputs from pet-volume pipeline
        "acq_label": "FDG",
        "use_pvc_data": False,
        "suvr_reference_region": "pons",
    }

    pipeline = StatisticsVolume(
        caps_directory=join(root, "in", "caps"),
        tsv_file=join(root, "in", "group-UnitTest_covariates.tsv"),
        parameters=parameters,
    )
    pipeline.build()


def test_instantiate_StatisticsVolumeCorrection():
    from os.path import abspath, dirname, join

    from clinica.pipelines.statistics_volume_correction.statistics_volume_correction_pipeline import (
        StatisticsVolumeCorrection,
    )

    root = dirname(abspath(join(abspath(__file__), pardir)))
    root = join(root, "data", "StatisticsVolumeCorrection")

    # Instantiate pipeline and run()
    parameters = {
        "t_map": "group-UnitTest_AD-lt-CN_measure-fdg_fwhm-8_TStatistics.nii",
        "height_threshold": 3.2422,
        "FWEp": 4.928,
        "FDRp": 4.693,
        "FWEc": 206987,
        "FDRc": 206987,
        "n_cuts": 15,
    }
    pipeline = StatisticsVolumeCorrection(
        caps_directory=join(root, "in", "caps"), parameters=parameters
    )
    pipeline.build()
