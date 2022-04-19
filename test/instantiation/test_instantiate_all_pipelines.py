# coding: utf8

import warnings
from os import fspath
from pathlib import Path

# Small unit tests for all pipelines
##
# test if instantiation and building of workflows is working


warnings.filterwarnings("ignore")


def test_instantiate_T1FreeSurferCrossSectional(cmdopt):

    from clinica.pipelines.t1_freesurfer.t1_freesurfer_pipeline import T1FreeSurfer

    input_dir = Path(cmdopt["input"])
    root = input_dir / "T1FreeSurfer"

    parameters = {
        "recon_all_args": "-qcache",
        "skip_question": True,
    }
    pipeline = T1FreeSurfer(
        bids_directory=fspath(root / "in" / "bids"),
        caps_directory=fspath(root / "in" / "caps"),
        tsv_file=fspath(root / "in" / "subjects.tsv"),
        parameters=parameters,
    )
    pipeline.build()


def test_instantiate_T1VolumeTissueSegmentation(cmdopt):

    from clinica.pipelines.t1_volume_tissue_segmentation.t1_volume_tissue_segmentation_pipeline import (
        T1VolumeTissueSegmentation,
    )

    input_dir = Path(cmdopt["input"])
    root = input_dir / "T1VolumeTissueSegmentation"

    parameters = {"skip_question": True}
    pipeline = T1VolumeTissueSegmentation(
        bids_directory=fspath(root / "in" / "bids"),
        caps_directory=fspath(root / "in" / "caps"),
        tsv_file=fspath(root / "in" / "subjects.tsv"),
        parameters=parameters,
    )
    pipeline.build()


def test_instantiate_T1VolumeCreateDartel(cmdopt):

    from clinica.pipelines.t1_volume_create_dartel.t1_volume_create_dartel_pipeline import (
        T1VolumeCreateDartel,
    )

    input_dir = Path(cmdopt["input"])
    root = input_dir / "T1VolumeCreateDartel"

    parameters = {"group_label": "UnitTest"}
    pipeline = T1VolumeCreateDartel(
        bids_directory=fspath(root / "in" / "bids"),
        caps_directory=fspath(root / "in" / "caps"),
        tsv_file=fspath(root / "in" / "subjects.tsv"),
        parameters=parameters,
    )
    pipeline.build()


def test_instantiate_T1VolumeDartel2MNI(cmdopt):

    from clinica.pipelines.t1_volume_dartel2mni.t1_volume_dartel2mni_pipeline import (
        T1VolumeDartel2MNI,
    )

    input_dir = Path(cmdopt["input"])
    root = input_dir / "T1VolumeDartel2MNI"

    parameters = {"group_label": "UnitTest"}
    pipeline = T1VolumeDartel2MNI(
        bids_directory=fspath(root / "in" / "bids"),
        caps_directory=fspath(root / "in" / "caps"),
        tsv_file=fspath(root / "in" / "subjects.tsv"),
        parameters=parameters,
    )
    pipeline.build()


def test_instantiate_T1VolumeRegisterDartel(cmdopt):

    from clinica.pipelines.t1_volume_register_dartel.t1_volume_register_dartel_pipeline import (
        T1VolumeRegisterDartel,
    )

    input_dir = Path(cmdopt["input"])
    root = input_dir / "T1VolumeRegisterDartel"

    parameters = {"group_label": "UnitTest"}
    pipeline = T1VolumeRegisterDartel(
        bids_directory=fspath(root / "in" / "bids"),
        caps_directory=fspath(root / "in" / "caps"),
        tsv_file=fspath(root / "in" / "subjects.tsv"),
        parameters=parameters,
    )
    pipeline.build()


def test_instantiate_T1VolumeParcellation(cmdopt):

    from clinica.pipelines.t1_volume_parcellation.t1_volume_parcellation_pipeline import (
        T1VolumeParcellation,
    )

    input_dir = Path(cmdopt["input"])
    root = input_dir / "T1VolumeParcellation"

    parameters = {"group_label": "UnitTest"}
    pipeline = T1VolumeParcellation(
        caps_directory=fspath(root / "in" / "caps"),
        tsv_file=fspath(root / "in" / "subjects.tsv"),
        parameters=parameters,
    )
    pipeline.build()


def test_instantiate_DWIPreprocessingUsingT1(cmdopt):

    from clinica.pipelines.dwi_preprocessing_using_t1.dwi_preprocessing_using_t1_pipeline import (
        DwiPreprocessingUsingT1,
    )

    input_dir = Path(cmdopt["input"])
    root = input_dir / "DWIPreprocessingUsingT1"

    parameters = {
        "initrand": False,
        "low_bval": 5,
        "use_cuda": False,
    }
    pipeline = DwiPreprocessingUsingT1(
        bids_directory=fspath(root / "in" / "bids"),
        caps_directory=fspath(root / "in" / "caps"),
        tsv_file=fspath(root / "in" / "subjects.tsv"),
        parameters=parameters,
    )
    pipeline.build()


def test_instantiate_DWIPreprocessingUsingPhaseDiffFieldmap(cmdopt):

    from clinica.pipelines.dwi_preprocessing_using_fmap.dwi_preprocessing_using_phasediff_fmap_pipeline import (
        DwiPreprocessingUsingPhaseDiffFMap,
    )

    input_dir = Path(cmdopt["input"])
    root = input_dir / "DWIPreprocessingUsingPhaseDiffFieldmap"

    parameters = {
        "initrand": False,
        "low_bval": 5,
        "use_cuda": False,
    }
    pipeline = DwiPreprocessingUsingPhaseDiffFMap(
        bids_directory=fspath(root / "in" / "bids"),
        caps_directory=fspath(root / "in" / "caps"),
        tsv_file=fspath(root / "in" / "subjects.tsv"),
        parameters=parameters,
    )
    pipeline.build()


def test_instantiate_DWIDTI(cmdopt):

    from clinica.pipelines.dwi_dti.dwi_dti_pipeline import DwiDti

    input_dir = Path(cmdopt["input"])
    root = input_dir / "DWIDTI"

    pipeline = DwiDti(
        caps_directory=fspath(root / "in" / "caps"),
        tsv_file=fspath(root / "in" / "subjects.tsv"),
    )
    pipeline.build()


def test_instantiate_DWIConnectome(cmdopt):

    from clinica.pipelines.dwi_connectome.dwi_connectome_pipeline import DwiConnectome

    input_dir = Path(cmdopt["input"])
    root = input_dir / "DWIConnectome"

    parameters = {"n_tracks": 1000}
    pipeline = DwiConnectome(
        caps_directory=fspath(root / "in" / "caps"),
        tsv_file=fspath(root / "in" / "subjects.tsv"),
        parameters=parameters,
    )
    pipeline.build()


def test_instantiate_PETVolume(cmdopt):

    from clinica.pipelines.pet_volume.pet_volume_pipeline import PETVolume
    from clinica.utils.pet import Tracer

    input_dir = Path(cmdopt["input"])
    root = input_dir / "PETVolume"

    parameters = {
        "group_label": "UnitTest",
        "acq_label": Tracer.FDG,
        "suvr_reference_region": "pons",
        "skip_question": True,
    }
    pipeline = PETVolume(
        bids_directory=fspath(root / "in" / "bids"),
        caps_directory=fspath(root / "in" / "caps"),
        tsv_file=fspath(root / "in" / "subjects.tsv"),
        parameters=parameters,
    )
    pipeline.build()


def test_instantiate_PETLinear(cmdopt):

    from clinica.pipelines.pet_linear.pet_linear_pipeline import PETLinear
    from clinica.utils.pet import Tracer

    input_dir = Path(cmdopt["input"])
    root = input_dir / "PETLinear"

    parameters = {
        "acq_label": Tracer.FDG,
        "suvr_reference_region": "cerebellumPons2",
        "skip_question": True,
    }
    pipeline = PETLinear(
        bids_directory=fspath(root / "in" / "bids"),
        caps_directory=fspath(root / "in" / "caps"),
        tsv_file=fspath(root / "in" / "subjects.tsv"),
        parameters=parameters,
    )
    pipeline.build


def test_instantiate_StatisticsSurface(cmdopt):

    from clinica.pipelines.statistics_surface.statistics_surface_pipeline import (
        StatisticsSurface,
    )

    input_dir = Path(cmdopt["input"])
    root = input_dir / "StatisticsSurface"

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
        caps_directory=fspath(root / "in" / "caps"),
        tsv_file=fspath(root / "in" / "subjects.tsv"),
        parameters=parameters,
    )

    pipeline.build()


def test_instantiate_PETSurfaceCrossSectional(cmdopt):

    from clinica.pipelines.pet_surface.pet_surface_pipeline import PetSurface
    from clinica.utils.pet import Tracer

    input_dir = Path(cmdopt["input"])
    root = input_dir / "PETSurface"

    parameters = {
        "acq_label": Tracer.FDG,
        "suvr_reference_region": "pons",
        "pvc_psf_tsv": fspath(root / "in" / "subjects.tsv"),
        "longitudinal": False,
        "skip_question": True,
    }
    pipeline = PetSurface(
        bids_directory=fspath(root / "in" / "bids"),
        caps_directory=fspath(root / "in" / "caps"),
        tsv_file=fspath(root / "in" / "subjects.tsv"),
        parameters=parameters,
    )
    pipeline.build()


# def test_instantiate_PETSurfaceLongitudinal():
#     from os.path import dirname, join, abspath
#     from clinica.pipelines.pet_surface.pet_surface_pipeline import PetSurface
#     from clinica.utils.pet import Tracer
#
#     root = dirname(abspath(join(abspath(__file__), pardir)))
#     root = join(root, 'data', 'PETSurfaceLongitudinal')
#     parameters = {
#         'acq_label': Tracer.FDG,
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


def test_instantiate_WorlflowsML(cmdopt):
    from os.path import exists

    from clinica.pipelines.machine_learning.input import (
        CAPSRegionBasedInput,
        CAPSVertexBasedInput,
        CAPSVoxelBasedInput,
    )
    from clinica.utils.pet import Tracer

    input_dir = Path(cmdopt["input"])
    root = input_dir / "WorkflowsML"

    caps_dir = fspath(root / "in" / "caps")
    tsv = fspath(root / "in" / "subjects.tsv")
    diagnoses_tsv = fspath(root / "in" / "diagnosis.tsv")
    group_label = "allADNIdartel"
    image_type = ["T1w", "PET"]
    atlases = ["AAL2", "Neuromorphometrics", "AICHA", "LPBA40", "Hammers"]
    possible_fwhm = [0, 5, 10, 15, 20, 25]
    tracer = Tracer.FDG

    voxel_input = [
        CAPSVoxelBasedInput(
            {
                "caps_directory": caps_dir,
                "subjects_visits_tsv": tsv,
                "diagnoses_tsv": diagnoses_tsv,
                "group_label": group_label,
                "image_type": im,
                "fwhm": 8,
                "acq_label": tracer,
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
                "acq_label": tracer,
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
                "acq_label": tracer,
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


def test_instantiate_SpatialSVM(cmdopt):

    from clinica.pipelines.machine_learning_spatial_svm.spatial_svm_pipeline import (
        SpatialSVM,
    )

    input_dir = Path(cmdopt["input"])
    root = input_dir / "SpatialSVM"

    parameters = {"group_label": "ADNIbl", "orig_input_data": "t1-volume"}
    pipeline = SpatialSVM(
        caps_directory=fspath(root / "in" / "caps"),
        tsv_file=fspath(root / "in" / "subjects.tsv"),
        parameters=parameters,
    )
    pipeline.build()


def test_instantiate_T1FreeSurferTemplate(cmdopt):

    from clinica.pipelines.t1_freesurfer_longitudinal.t1_freesurfer_template_pipeline import (
        T1FreeSurferTemplate,
    )

    input_dir = Path(cmdopt["input"])
    root = input_dir / "T1FreeSurferTemplate"

    pipeline = T1FreeSurferTemplate(
        caps_directory=fspath(root / "in" / "caps"),
        tsv_file=fspath(root / "in" / "subjects.tsv"),
    )
    pipeline.build()


def test_instantiate_T1FreeSurferLongitudinalCorrection(cmdopt):

    from clinica.pipelines.t1_freesurfer_longitudinal.t1_freesurfer_longitudinal_correction_pipeline import (
        T1FreeSurferLongitudinalCorrection,
    )

    input_dir = Path(cmdopt["input"])
    root = input_dir / "T1FreeSurferLongitudinalCorrection"

    pipeline = T1FreeSurferLongitudinalCorrection(
        caps_directory=fspath(root / "in" / "caps"),
        tsv_file=fspath(root / "in" / "subjects.tsv"),
    )
    pipeline.build()


def test_instantiate_T1Linear(cmdopt):

    from clinica.pipelines.t1_linear.t1_linear_pipeline import T1Linear

    input_dir = Path(cmdopt["input"])
    root = input_dir / "T1Linear"

    parameters = {"uncropped_image": False}

    pipeline = T1Linear(
        bids_directory=fspath(root / "in" / "bids"),
        caps_directory=fspath(root / "in" / "caps"),
        tsv_file=fspath(root / "in" / "subjects.tsv"),
    )
    pipeline.build()


def test_instantiate_StatisticsVolume(cmdopt):

    from clinica.pipelines.statistics_volume.statistics_volume_pipeline import (
        StatisticsVolume,
    )
    from clinica.utils.pet import Tracer

    input_dir = Path(cmdopt["input"])
    root = input_dir / "StatisticsVolume"

    # Instantiate pipeline and run()
    parameters = {
        # Clinica compulsory parameters
        "group_label": "UnitTest",
        "orig_input_data_volume": "pet-volume",
        "contrast": "group",
        # Optional arguments for inputs from pet-volume pipeline
        "acq_label": Tracer.FDG,
        "use_pvc_data": False,
        "suvr_reference_region": "pons",
    }

    pipeline = StatisticsVolume(
        caps_directory=fspath(root / "in" / "caps"),
        tsv_file=fspath(root / "in" / "group-UnitTest_covariates.tsv"),
        parameters=parameters,
    )
    pipeline.build()


def test_instantiate_StatisticsVolumeCorrection(cmdopt):

    from clinica.pipelines.statistics_volume_correction.statistics_volume_correction_pipeline import (
        StatisticsVolumeCorrection,
    )

    input_dir = Path(cmdopt["input"])
    root = input_dir / "StatisticsVolumeCorrection"

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
        caps_directory=fspath(root / "in" / "caps"), parameters=parameters
    )
    pipeline.build()
