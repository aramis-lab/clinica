import shutil
import warnings
from os import fspath
from pathlib import Path
from test.nonregression.testing_tools import configure_paths

import pytest

from clinica.utils.pet import SUVRReferenceRegion, Tracer

warnings.filterwarnings("ignore")


def test_instantiate_t1_freesurfer_cross_sectional(cmdopt, tmp_path):
    from clinica.pipelines.anatomical.freesurfer.t1.pipeline import T1FreeSurfer

    base_dir = Path(cmdopt["input"])
    working_dir = Path(cmdopt["wd"])
    input_dir, tmp_dir, ref_dir = configure_paths(base_dir, tmp_path, "T1FreeSurfer")
    T1FreeSurfer(
        bids_directory=fspath(input_dir / "bids"),
        caps_directory=fspath(tmp_dir / "caps"),
        tsv_file=fspath(input_dir / "subjects.tsv"),
        base_dir=fspath(working_dir),
        parameters={
            "recon_all_args": "-qcache",
            "skip_question": True,
        },
    ).build()


def test_instantiate_spm_based_pipelines(cmdopt, tmp_path):
    """Run the instantiation tests for all pipelines using SPM.
    This avoids issues when running tests in parallel.
    """
    base_dir = Path(cmdopt["input"])
    working_dir = Path(cmdopt["wd"])

    run_tissue_segmentation(base_dir, tmp_path, working_dir)
    run_create_dartel(base_dir, tmp_path, working_dir)
    run_dartel_to_mni(base_dir, tmp_path, working_dir)
    run_register_dartel(base_dir, tmp_path, working_dir)
    run_pet_volume(base_dir, tmp_path, working_dir)


def run_tissue_segmentation(
    base_dir: Path, output_dir: Path, working_dir: Path
) -> None:
    from clinica.pipelines.t1_volume_tissue_segmentation.t1_volume_tissue_segmentation_pipeline import (
        T1VolumeTissueSegmentation,
    )

    input_dir, tmp_dir, ref_dir = configure_paths(
        base_dir, output_dir, "T1VolumeTissueSegmentation"
    )
    T1VolumeTissueSegmentation(
        bids_directory=fspath(input_dir / "bids"),
        caps_directory=fspath(tmp_dir / "caps"),
        tsv_file=fspath(input_dir / "subjects.tsv"),
        base_dir=fspath(working_dir),
        parameters={"skip_question": True},
    ).build()


def run_create_dartel(base_dir: Path, output_dir: Path, working_dir: Path) -> None:
    from clinica.pipelines.t1_volume_create_dartel.t1_volume_create_dartel_pipeline import (
        T1VolumeCreateDartel,
    )

    input_dir, tmp_dir, ref_dir = configure_paths(
        base_dir, output_dir, "T1VolumeCreateDartel"
    )
    # Copy the CAPS folder in temp folder in order to have writing privileges
    shutil.copytree(str(input_dir / "caps"), str(tmp_dir / "caps"))
    T1VolumeCreateDartel(
        bids_directory=fspath(input_dir / "bids"),
        caps_directory=fspath(tmp_dir / "caps"),
        tsv_file=fspath(input_dir / "subjects.tsv"),
        base_dir=fspath(working_dir),
        group_label="UnitTest",
    ).build()


def run_dartel_to_mni(base_dir: Path, output_dir: Path, working_dir: Path) -> None:
    from clinica.pipelines.t1_volume_dartel2mni.t1_volume_dartel2mni_pipeline import (
        T1VolumeDartel2MNI,
    )

    input_dir, tmp_dir, ref_dir = configure_paths(
        base_dir, output_dir, "T1VolumeDartel2MNI"
    )
    # Copy the CAPS folder in temp folder in order to have writing privileges
    shutil.copytree(str(input_dir / "caps"), str(tmp_dir / "caps"))
    T1VolumeDartel2MNI(
        bids_directory=fspath(input_dir / "bids"),
        caps_directory=fspath(tmp_dir / "caps"),
        tsv_file=fspath(input_dir / "subjects.tsv"),
        base_dir=fspath(working_dir),
        group_label="UnitTest",
    ).build()


def run_register_dartel(base_dir: Path, output_dir: Path, working_dir: Path) -> None:
    from clinica.pipelines.t1_volume_register_dartel.t1_volume_register_dartel_pipeline import (
        T1VolumeRegisterDartel,
    )

    input_dir, tmp_dir, ref_dir = configure_paths(
        base_dir, output_dir, "T1VolumeRegisterDartel"
    )
    # Copy the CAPS folder in temp folder in order to have writing privileges
    shutil.copytree(str(input_dir / "caps"), str(tmp_dir / "caps"))
    T1VolumeRegisterDartel(
        bids_directory=fspath(input_dir / "bids"),
        caps_directory=fspath(tmp_dir / "caps"),
        tsv_file=fspath(input_dir / "subjects.tsv"),
        base_dir=fspath(working_dir),
        group_label="UnitTest",
    ).build()


def test_instantiate_t1_volume_parcellation(cmdopt, tmp_path):
    from clinica.pipelines.t1_volume_parcellation.t1_volume_parcellation_pipeline import (
        T1VolumeParcellation,
    )

    base_dir = Path(cmdopt["input"])
    working_dir = Path(cmdopt["wd"])
    input_dir, tmp_dir, ref_dir = configure_paths(
        base_dir, tmp_path, "T1VolumeParcellation"
    )
    # Copy the CAPS folder in temp folder in order to have writing privileges
    shutil.copytree(str(input_dir / "caps"), str(tmp_dir / "caps"))
    T1VolumeParcellation(
        caps_directory=fspath(tmp_dir / "caps"),
        tsv_file=fspath(input_dir / "subjects.tsv"),
        base_dir=fspath(working_dir),
        group_label="UnitTest",
    ).build()


def test_instantiate_dwi_preprocessing_using_t1(cmdopt, tmp_path):
    from clinica.pipelines.dwi.preprocessing.t1.pipeline import (
        DwiPreprocessingUsingT1,
    )

    base_dir = Path(cmdopt["input"])
    working_dir = Path(cmdopt["wd"])
    input_dir, tmp_dir, ref_dir = configure_paths(
        base_dir, tmp_path, "DWIPreprocessingUsingT1"
    )
    DwiPreprocessingUsingT1(
        bids_directory=fspath(input_dir / "bids"),
        caps_directory=fspath(tmp_dir / "caps"),
        tsv_file=fspath(input_dir / "subjects.tsv"),
        base_dir=fspath(working_dir),
        parameters={
            "initrand": False,
            "low_bval": 5,
            "use_cuda": False,
            "delete_cache": True,
            "random_seed": 42,
        },
    ).build()


def test_instantiate_dwi_preprocessing_using_phase_diff_field_map(cmdopt, tmp_path):
    from clinica.pipelines.dwi.preprocessing.fmap.pipeline import (
        DwiPreprocessingUsingPhaseDiffFMap,
    )

    base_dir = Path(cmdopt["input"])
    working_dir = Path(cmdopt["wd"])
    input_dir, tmp_dir, ref_dir = configure_paths(
        base_dir, tmp_path, "DWIPreprocessingUsingPhaseDiffFieldmap"
    )
    DwiPreprocessingUsingPhaseDiffFMap(
        bids_directory=fspath(input_dir / "bids"),
        caps_directory=fspath(tmp_dir / "caps"),
        tsv_file=fspath(input_dir / "subjects.tsv"),
        base_dir=fspath(working_dir),
        parameters={
            "initrand": False,
            "low_bval": 5,
            "use_cuda": False,
        },
    ).build()


def test_instantiate_dwi_dti(cmdopt, tmp_path):
    from clinica.pipelines.dwi.dti.pipeline import DwiDti

    base_dir = Path(cmdopt["input"])
    working_dir = Path(cmdopt["wd"])
    input_dir, tmp_dir, ref_dir = configure_paths(base_dir, tmp_path, "DWIDTI")
    # Copy the CAPS folder in temp folder in order to have writing privileges
    shutil.copytree(str(input_dir / "caps"), str(tmp_dir / "caps"))
    DwiDti(
        caps_directory=fspath(tmp_dir / "caps"),
        tsv_file=fspath(input_dir / "subjects.tsv"),
        base_dir=fspath(working_dir),
    ).build()


def test_instantiate_dwi_connectome(cmdopt, tmp_path):
    from clinica.pipelines.dwi.connectome.pipeline import DwiConnectome

    base_dir = Path(cmdopt["input"])
    working_dir = Path(cmdopt["wd"])
    input_dir, tmp_dir, ref_dir = configure_paths(base_dir, tmp_path, "DWIConnectome")
    # Copy the CAPS folder in temp folder in order to have writing privileges
    shutil.copytree(str(input_dir / "caps"), str(tmp_dir / "caps"))
    DwiConnectome(
        caps_directory=fspath(tmp_dir / "caps"),
        tsv_file=fspath(input_dir / "subjects.tsv"),
        base_dir=fspath(working_dir),
        parameters={"n_tracks": 1000},
    ).build()


def run_pet_volume(base_dir: Path, output_dir: Path, working_dir: Path) -> None:
    from clinica.pipelines.pet.volume.pipeline import PETVolume

    input_dir, tmp_dir, ref_dir = configure_paths(base_dir, output_dir, "PETVolume")
    # Copy the CAPS folder in temp folder in order to have writing privileges
    shutil.copytree(str(input_dir / "caps"), str(tmp_dir / "caps"))
    PETVolume(
        bids_directory=fspath(input_dir / "bids"),
        caps_directory=fspath(tmp_dir / "caps"),
        tsv_file=fspath(input_dir / "subjects.tsv"),
        base_dir=fspath(working_dir),
        group_label="UnitTest",
        parameters={
            "acq_label": Tracer.FDG,
            "suvr_reference_region": SUVRReferenceRegion.PONS,
            "skip_question": True,
            "reconstruction_method": None,
        },
    ).build()


def test_instantiate_pet_linear(cmdopt, tmp_path):
    from clinica.pipelines.pet.linear.pipeline import PETLinear

    base_dir = Path(cmdopt["input"])
    working_dir = Path(cmdopt["wd"])
    input_dir, tmp_dir, ref_dir = configure_paths(base_dir, tmp_path, "PETLinear")
    # Copy the CAPS folder in temp folder in order to have writing privileges
    shutil.copytree(str(input_dir / "caps"), str(tmp_dir / "caps"))
    PETLinear(
        bids_directory=fspath(input_dir / "bids"),
        caps_directory=fspath(tmp_dir / "caps"),
        tsv_file=fspath(input_dir / "subjects.tsv"),
        base_dir=fspath(working_dir),
        parameters={
            "acq_label": Tracer.FDG,
            "suvr_reference_region": SUVRReferenceRegion.CEREBELLUM_PONS2,
            "skip_question": True,
            "reconstruction_method": None,
        },
    ).build()


def test_instantiate_statistics_surface(cmdopt, tmp_path):
    from clinica.pipelines.statistics_surface.pipeline import StatisticsSurface

    base_dir = Path(cmdopt["input"])
    working_dir = Path(cmdopt["wd"])
    input_dir, tmp_dir, ref_dir = configure_paths(
        base_dir, tmp_path, "StatisticsSurface"
    )
    # Copy the CAPS folder in temp folder in order to have writing privileges
    shutil.copytree(str(input_dir / "caps"), str(tmp_dir / "caps"))
    StatisticsSurface(
        caps_directory=fspath(tmp_dir / "caps"),
        tsv_file=fspath(input_dir / "subjects.tsv"),
        base_dir=fspath(working_dir),
        group_label="UnitTest",
        parameters={
            "orig_input_data": "t1-freesurfer",
            "glm_type": "group_comparison",
            "contrast": "group",
            "covariates": "age sex",
        },
    ).build()


def test_instantiate_pet_surface_cross_sectional(cmdopt, tmp_path):
    from clinica.pipelines.pet_surface.pet_surface_pipeline import PetSurface

    base_dir = Path(cmdopt["input"])
    working_dir = Path(cmdopt["wd"])
    input_dir, tmp_dir, ref_dir = configure_paths(base_dir, tmp_path, "PETSurface")
    # Copy the CAPS folder in temp folder in order to have writing privileges
    shutil.copytree(str(input_dir / "caps"), str(tmp_dir / "caps"))
    PetSurface(
        bids_directory=fspath(input_dir / "bids"),
        caps_directory=fspath(tmp_dir / "caps"),
        tsv_file=fspath(input_dir / "subjects.tsv"),
        base_dir=fspath(working_dir),
        parameters={
            "acq_label": Tracer.FDG,
            "suvr_reference_region": SUVRReferenceRegion.PONS,
            "pvc_psf_tsv": fspath(input_dir / "subjects.tsv"),
            "longitudinal": False,
            "skip_question": True,
            "reconstruction_method": None,
        },
    ).build()


@pytest.mark.skip(reason="Currently broken. Needs to be fixed...")
def test_instantiate_pet_surface_longitudinal(cmdopt):
    from clinica.pipelines.pet_surface.pet_surface_pipeline import PetSurface

    input_dir = Path(cmdopt["input"])
    root = input_dir / "PETSurfaceLongitudinal"
    parameters = {
        "acq_label": Tracer.FDG,
        "suvr_reference_region": SUVRReferenceRegion.PONS,
        "pvc_psf_tsv": fspath(root / "in" / "subjects.tsv"),
        "longitudinal": True,
        "reconstruction_method": None,
    }
    pipeline = PetSurface(
        bids_directory=fspath(root / "in" / "bids"),
        caps_directory=fspath(root / "in" / "caps"),
        tsv_file=fspath(root / "in" / "subjects.tsv"),
        parameters=parameters,
    )
    pipeline.build()


def test_instantiate_workflows_ml(cmdopt, tmp_path):
    from clinica.pipelines.machine_learning.input import (
        CAPSRegionBasedInput,
        CAPSVertexBasedInput,
        CAPSVoxelBasedInput,
    )

    base_dir = Path(cmdopt["input"])
    working_dir = Path(cmdopt["wd"])
    input_dir, tmp_dir, ref_dir = configure_paths(base_dir, tmp_path, "WorkflowsML")
    # Copy the CAPS folder in temp folder in order to have writing privileges
    shutil.copytree(str(input_dir / "caps"), str(tmp_dir / "caps"))
    caps_dir = fspath(tmp_dir / "caps")
    tsv = fspath(input_dir / "subjects.tsv")
    diagnoses_tsv = fspath(input_dir / "diagnosis.tsv")
    group_label = "allADNIdartel"
    image_type = ["T1w", "PET"]
    atlases = ["AAL2", "Neuromorphometrics", "AICHA", "LPBA40", "Hammers"]
    possible_fwhm = [0, 5, 10, 15, 20, 25]
    tracer = Tracer.FDG
    region = SUVRReferenceRegion.PONS
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
                "suvr_reference_region": region,
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
                "suvr_reference_region": region,
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
                "suvr_reference_region": region,
            }
        )
        for fwhm in possible_fwhm
    ]
    # Check that each file exists
    for inputs in voxel_input + region_input + vertex_input:
        for file in inputs.get_images():
            if isinstance(file, str):
                assert Path(file).exists()
            elif isinstance(file, list) and len(file) == 2:
                assert all([Path(p).exists() for p in file])
            else:
                raise ValueError("An error occurred...")


def test_instantiate_spatial_svm(cmdopt, tmp_path):
    from clinica.pipelines.machine_learning_spatial_svm.spatial_svm_pipeline import (
        SpatialSVM,
    )

    base_dir = Path(cmdopt["input"])
    working_dir = Path(cmdopt["wd"])
    input_dir, tmp_dir, ref_dir = configure_paths(base_dir, tmp_path, "SpatialSVM")
    # Copy the CAPS folder in temp folder in order to have writing privileges
    shutil.copytree(str(input_dir / "caps"), str(tmp_dir / "caps"))
    SpatialSVM(
        caps_directory=fspath(tmp_dir / "caps"),
        tsv_file=fspath(input_dir / "subjects.tsv"),
        base_dir=fspath(working_dir),
        group_label="ADNIbl",
        parameters={"orig_input_data_ml": "t1-volume"},
    ).build()


def test_instantiate_t1_freesurfer_template(cmdopt, tmp_path):
    from clinica.pipelines.anatomical.freesurfer.longitudinal.template.pipeline import (
        T1FreeSurferTemplate,
    )

    base_dir = Path(cmdopt["input"])
    working_dir = Path(cmdopt["wd"])
    input_dir, tmp_dir, ref_dir = configure_paths(
        base_dir, tmp_path, "T1FreeSurferTemplate"
    )
    # Copy the CAPS folder in temp folder in order to have writing privileges
    shutil.copytree(str(input_dir / "caps"), str(tmp_dir / "caps"))
    T1FreeSurferTemplate(
        caps_directory=fspath(tmp_dir / "caps"),
        tsv_file=fspath(input_dir / "subjects.tsv"),
        base_dir=fspath(working_dir),
    ).build()


def test_instantiate_t1_freesurfer_longitudinal_correction(cmdopt, tmp_path):
    from clinica.pipelines.anatomical.freesurfer.longitudinal.correction.pipeline import (
        T1FreeSurferLongitudinalCorrection,
    )

    base_dir = Path(cmdopt["input"])
    working_dir = Path(cmdopt["wd"])
    input_dir, tmp_dir, ref_dir = configure_paths(
        base_dir, tmp_path, "T1FreeSurferLongitudinalCorrection"
    )
    # Copy the CAPS folder in temp folder in order to have writing privileges
    shutil.copytree(str(input_dir / "caps"), str(tmp_dir / "caps"))
    T1FreeSurferLongitudinalCorrection(
        caps_directory=fspath(tmp_dir / "caps"),
        tsv_file=fspath(input_dir / "subjects.tsv"),
        base_dir=fspath(working_dir),
    ).build()


def test_instantiate_t1_linear(cmdopt, tmp_path):
    from clinica.pipelines.t1_linear.anat_linear_pipeline import AnatLinear

    base_dir = Path(cmdopt["input"])
    working_dir = Path(cmdopt["wd"])
    input_dir, tmp_dir, ref_dir = configure_paths(base_dir, tmp_path, "T1Linear")
    AnatLinear(
        bids_directory=fspath(input_dir / "bids"),
        caps_directory=fspath(tmp_dir / "caps"),
        tsv_file=fspath(input_dir / "subjects.tsv"),
        base_dir=fspath(working_dir),
        name="t1-linear",
        parameters={"uncropped_image": False},
    ).build()


def test_instantiate_statistics_volume(cmdopt, tmp_path):
    from clinica.pipelines.statistics_volume.statistics_volume_pipeline import (
        StatisticsVolume,
    )

    base_dir = Path(cmdopt["input"])
    working_dir = Path(cmdopt["wd"])
    input_dir, tmp_dir, ref_dir = configure_paths(
        base_dir, tmp_path, "StatisticsVolume"
    )
    # Copy the CAPS folder in temp folder in order to have writing privileges
    shutil.copytree(str(input_dir / "caps"), str(tmp_dir / "caps"))

    StatisticsVolume(
        caps_directory=fspath(tmp_dir / "caps"),
        tsv_file=fspath(input_dir / "group-UnitTest_covariates.tsv"),
        base_dir=fspath(working_dir),
        group_label="UnitTest",
        parameters={
            "orig_input_data_volume": "pet-volume",
            "contrast": "group",
            "acq_label": Tracer.FDG,
            "use_pvc_data": False,
            "suvr_reference_region": SUVRReferenceRegion.PONS,
        },
    ).build()


def test_instantiate_statistics_volume_correction(cmdopt, tmp_path):
    from clinica.pipelines.statistics_volume_correction.statistics_volume_correction_pipeline import (
        StatisticsVolumeCorrection,
    )

    base_dir = Path(cmdopt["input"])
    working_dir = Path(cmdopt["wd"])
    input_dir, tmp_dir, ref_dir = configure_paths(
        base_dir, tmp_path, "StatisticsVolumeCorrection"
    )
    # Copy the CAPS folder in temp folder in order to have writing privileges
    shutil.copytree(str(input_dir / "caps"), str(tmp_dir / "caps"))
    StatisticsVolumeCorrection(
        caps_directory=fspath(tmp_dir / "caps"),
        base_dir=fspath(working_dir),
        parameters={
            "t_map": "group-UnitTest_AD-lt-CN_measure-fdg_fwhm-8_TStatistics.nii",
            "height_threshold": 3.2422,
            "FWEp": 4.928,
            "FDRp": 4.693,
            "FWEc": 206987,
            "FDRc": 206987,
            "n_cuts": 15,
        },
    ).build()
