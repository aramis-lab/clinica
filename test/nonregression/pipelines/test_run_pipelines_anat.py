# coding: utf8

"""
This file contains a set of functional tests designed to check the correct execution of the pipeline and the
different functions available in Clinica
"""

import warnings
from os import fspath
from pathlib import Path
from test.nonregression.testing_tools import *

import pytest
import functools

# Determine location for working_directory
warnings.filterwarnings("ignore")


def configure_paths(base_dir, tmp_path, name):
    input_dir = base_dir / name / "in"
    ref_dir = base_dir / name / "ref"
    tmp_out_dir = tmp_path / name / "out"
    tmp_out_dir.mkdir(parents=True)
    
    return input_dir, tmp_out_dir, ref_dir


@pytest.mark.fast
def test_t1_linear(cmdopt, tmp_path):
    base_dir = Path(cmdopt["input"])
    working_dir = Path(cmdopt["wd"])
    input_dir, tmp_dir, ref_dir = configure_paths(base_dir, tmp_path, "T1Linear")
    run_T1Linear(input_dir, tmp_dir, ref_dir, working_dir)


@pytest.mark.fast
def test_flair_linear(cmdopt, tmp_path):
    base_dir = Path(cmdopt["input"])
    working_dir = Path(cmdopt["wd"])
    input_dir, tmp_dir, ref_dir = configure_paths(base_dir, tmp_path, "FlairLinear")
    run_FlairLinear(input_dir, tmp_dir, ref_dir, working_dir)


@pytest.mark.slow
def test_t1_freesurfer_cross_sectional(cmdopt, tmp_path):
    base_dir = Path(cmdopt["input"])
    working_dir = Path(cmdopt["wd"])
    input_dir, tmp_dir, ref_dir = configure_paths(base_dir, tmp_path, "T1FreeSurfer")
    run_T1FreeSurferCrossSectional(input_dir, tmp_dir, ref_dir, working_dir)


def test_t1_volume_tissue_segmentation(cmdopt, tmp_path):
    base_dir = Path(cmdopt["input"])
    working_dir = Path(cmdopt["wd"])
    input_dir, tmp_dir, ref_dir = configure_paths(base_dir, tmp_path, "T1VolumeTissueSegmentation")
    run_T1VolumeTissueSegmentation(input_dir, tmp_dir, ref_dir, working_dir)


def test_t1_volume_create_dartel(cmdopt, tmp_path):
    base_dir = Path(cmdopt["input"])
    working_dir = Path(cmdopt["wd"])
    input_dir, tmp_dir, ref_dir = configure_paths(base_dir, tmp_path, "T1VolumeCreateDartel")
    run_T1VolumeCreateDartel(input_dir, tmp_dir, ref_dir, working_dir)


def test_t1_volume_dartel2mni(cmdopt, tmp_path):
    base_dir = Path(cmdopt["input"])
    working_dir = Path(cmdopt["wd"])
    input_dir, tmp_dir, ref_dir = configure_paths(base_dir, tmp_path, "T1VolumeDartel2MNI")
    run_T1VolumeDartel2MNI(input_dir, tmp_dir, ref_dir, working_dir)


def test_t1_volume_register_dartel(cmdopt, tmp_path):
    base_dir = Path(cmdopt["input"])
    working_dir = Path(cmdopt["wd"])
    input_dir, tmp_dir, ref_dir = configure_paths(base_dir, tmp_path, "T1VolumeRegisterDartel")
    run_T1VolumeRegisterDartel(input_dir, tmp_dir, ref_dir, working_dir)


def test_t1_volume_parcellation(cmdopt, tmp_path):
    base_dir = Path(cmdopt["input"])
    working_dir = Path(cmdopt["wd"])
    input_dir, tmp_dir, ref_dir = configure_paths(base_dir, tmp_path, "T1VolumeParcellation")
    run_T1VolumeParcellation(input_dir, tmp_dir, ref_dir, working_dir)


@pytest.mark.slow
def test_t1_freesurfer_template(cmdopt, tmp_path):
    base_dir = Path(cmdopt["input"])
    working_dir = Path(cmdopt["wd"])
    input_dir, tmp_dir, ref_dir = configure_paths(base_dir, tmp_path, "T1FreeSurferTemplate")
    run_T1FreeSurferTemplate(input_dir, tmp_dir, ref_dir, working_dir)


@pytest.mark.slow
def test_t1_freesurfer_longitudinal_correction(cmdopt, tmp_path):
    base_dir = Path(cmdopt["input"])
    working_dir = Path(cmdopt["wd"])
    input_dir, tmp_dir, ref_dir = configure_paths(base_dir, tmp_path, "T1FreeSurferLongitudinalCorrection")
    run_T1FreeSurferLongitudinalCorrection(input_dir, tmp_dir, ref_dir, working_dir)


def run_T1FreeSurferCrossSectional(
    input_dir: Path, output_dir: Path, ref_dir: Path, working_dir: Path
) -> None:
    # Data for this functional test comes from https://openneuro.org/datasets/ds000204

    from clinica.pipelines.t1_freesurfer.t1_freesurfer_pipeline import T1FreeSurfer

    parameters = {"recon_all_args": "-qcache", "skip_question": False}

    pipeline = T1FreeSurfer(
        bids_directory=fspath(input_dir / "bids"),
        caps_directory=fspath(output_dir / "caps"),
        tsv_file=fspath(input_dir / "subjects.tsv"),
        parameters=parameters,
        base_dir=fspath(working_dir),
    )
    pipeline.run(bypass_check=True)

    # We only check that folders are the same meaning that FreeSurfer finished without error
    # surf/ folder is ignored because it contains sym links that makes hard to check with ref data
    # (sym links of ref data are ignored after rsync on CI machines)
    def path_to_caps_fs(part_id: str, sess_id: str) -> Path:

        output_folder = (
            Path("caps")
            / "subjects"
            / part_id
            / sess_id
            / "t1"
            / "freesurfer_cross_sectional"
        )
        return output_folder

    folder1 = path_to_caps_fs("sub-01", "ses-2011")
    compare_folders(
        output_dir / folder1 / "regional_measures",
        ref_dir / folder1 / "regional_measures",
        output_dir,
    )
    compare_folders(
        output_dir / folder1 / "sub-01_ses-2011" / "label",
        ref_dir / folder1 / "sub-01_ses-2011" / "label",
        output_dir,
    )
    compare_folders(
        output_dir / folder1 / "sub-01_ses-2011" / "mri",
        ref_dir / folder1 / "sub-01_ses-2011" / "mri",
        output_dir,
    )
    compare_folders(
        output_dir / folder1 / "sub-01_ses-2011" / "stats",
        ref_dir / folder1 / "sub-01_ses-2011" / "stats",
        output_dir,
    )


def run_T1VolumeTissueSegmentation(
    input_dir: Path, output_dir: Path, ref_dir: Path, working_dir: Path
) -> None:

    from clinica.pipelines.t1_volume_tissue_segmentation.t1_volume_tissue_segmentation_pipeline import (
        T1VolumeTissueSegmentation,
    )

    parameters = {"skip_question": False}
    pipeline = T1VolumeTissueSegmentation(
        bids_directory=fspath(input_dir / "bids"),
        caps_directory=fspath(output_dir / "caps"),
        tsv_file=fspath(input_dir / "subjects.tsv"),
        parameters=parameters,
        base_dir=fspath(working_dir),
    )
    pipeline.build()
    pipeline.run(bypass_check=True)

    out_file = fspath(
        output_dir
        / "caps"
        / "subjects"
        / "sub-ADNI011S4105"
        / "ses-M0000"
        / "t1"
        / "spm"
        / "segmentation"
        / "dartel_input"
        / "sub-ADNI011S4105_ses-M000_T1w_segm-graymatter_dartelinput.nii.gz",
    )

    ref_file = fspath(
        ref_dir
        / "caps"
        / "subjects"
        / "sub-ADNI011S4105"
        / "ses-M000"
        / "t1"
        / "spm"
        / "segmentation"
        / "dartel_input"
        / "sub-ADNI011S4105_ses-M000_T1w_segm-graymatter_dartelinput.nii.gz",
    )

    assert likeliness_measure(out_file, ref_file, (1e-1, 0.02), (0.4, 0.01))


def run_T1VolumeCreateDartel(
    input_dir: Path, output_dir: Path, ref_dir: Path, working_dir: Path
) -> None:
    import shutil

    from clinica.pipelines.t1_volume_create_dartel.t1_volume_create_dartel_pipeline import (
        T1VolumeCreateDartel,
    )

    # Copy necessary data from in to out
    shutil.copytree(input_dir / "caps", output_dir / "caps", copy_function=shutil.copy)

    parameters = {"group_label": "UnitTest"}
    # Instantiate pipeline
    pipeline = T1VolumeCreateDartel(
        bids_directory=fspath(input_dir / "bids"),
        caps_directory=fspath(output_dir / "caps"),
        tsv_file=fspath(input_dir / "subjects.tsv"),
        base_dir=fspath(working_dir),
        parameters=parameters,
    )
    pipeline.build()
    pipeline.run(plugin="MultiProc", plugin_args={"n_procs": 4}, bypass_check=True)

    # Check output vs ref
    out_data_template = fspath(
        output_dir / "caps/groups/group-UnitTest/t1/group-UnitTest_template.nii.gz"
    )
    ref_data_template = fspath(ref_dir / "group-UnitTest_template.nii.gz")
    assert likeliness_measure(
        out_data_template, ref_data_template, (1e-3, 0.1), (1e-2, 0.1)
    )

    subjects = [
        "sub-ADNI011S4105",
        "sub-ADNI023S4020",
        "sub-ADNI035S4082",
        "sub-ADNI128S4832",
    ]
    out_data_forward_def = [
        fspath(
            output_dir
            / "caps"
            / "subjects"
            / sub
            / "ses-M000"
            / "t1"
            / "spm"
            / "dartel"
            / "group-UnitTest"
            / (
                sub
                + "_ses-M000_T1w_target-UnitTest_transformation-forward_deformation.nii.gz"
            )
        )
        for sub in subjects
    ]
    ref_data_forward_def = [
        fspath(
            ref_dir
            / (
                sub
                + "_ses-M000_T1w_target-UnitTest_transformation-forward_deformation.nii.gz"
            )
        )
        for sub in subjects
    ]

    for i in range(len(out_data_forward_def)):
        assert likeliness_measure(
            out_data_forward_def[i], ref_data_forward_def[i], (1e-3, 0.25), (1e-2, 0.1)
        )


def run_T1VolumeDartel2MNI(
    input_dir: Path, output_dir: Path, ref_dir: Path, working_dir: Path
) -> None:
    import shutil

    from clinica.pipelines.t1_volume_dartel2mni.t1_volume_dartel2mni_pipeline import (
        T1VolumeDartel2MNI,
    )

    # Copy necessary data from in to out
    shutil.copytree(input_dir / "caps", output_dir / "caps", copy_function=shutil.copy)

    parameters = {"group_label": "UnitTest"}
    # Instantiate pipeline and run()
    pipeline = T1VolumeDartel2MNI(
        bids_directory=fspath(input_dir / "bids"),
        caps_directory=fspath(output_dir / "caps"),
        tsv_file=fspath(input_dir / "subjects.tsv"),
        base_dir=fspath(working_dir),
        parameters=parameters,
    )
    pipeline.build()
    pipeline.run(plugin="MultiProc", plugin_args={"n_procs": 4}, bypass_check=True)

    # Check output vs ref
    subjects = [
        "sub-ADNI011S4105",
        "sub-ADNI023S4020",
        "sub-ADNI035S4082",
        "sub-ADNI128S4832",
    ]
    out_data_GM_MNI = [
        fspath(
            output_dir
            / "caps"
            / "subjects"
            / sub
            / "ses-M000"
            / "t1"
            / "spm"
            / "dartel"
            / "group-UnitTest"
            / (
                sub
                + "_ses-M000_T1w_segm-graymatter_space-Ixi549Space_modulated-on_fwhm-8mm_probability.nii.gz"
            )
        )
        for sub in subjects
    ]
    ref_data_GM_MNI = [
        fspath(
            ref_dir
            / (
                sub
                + "_ses-M000_T1w_segm-graymatter_space-Ixi549Space_modulated-on_fwhm-8mm_probability.nii.gz"
            )
        )
        for sub in subjects
    ]
    for i in range(len(out_data_GM_MNI)):
        assert likeliness_measure(
            out_data_GM_MNI[i], ref_data_GM_MNI[i], (1e-4, 0.15), (1, 0.02)
        )


def run_T1VolumeRegisterDartel(
    input_dir: Path, output_dir: Path, ref_dir: Path, working_dir: Path
) -> None:
    import shutil
    from os.path import abspath, dirname, join

    from clinica.pipelines.t1_volume_register_dartel.t1_volume_register_dartel_pipeline import (
        T1VolumeRegisterDartel,
    )

    # Copy necessary data from in to out
    shutil.copytree(input_dir / "caps", output_dir / "caps", copy_function=shutil.copy)

    # Instantiate and run pipeline
    parameters = {"group_label": "UnitTest"}
    pipeline = T1VolumeRegisterDartel(
        bids_directory=fspath(input_dir / "bids"),
        caps_directory=fspath(output_dir / "caps"),
        tsv_file=fspath(input_dir / "subjects.tsv"),
        base_dir=fspath(working_dir),
        parameters=parameters,
    )
    pipeline.build()
    pipeline.run(plugin="MultiProc", plugin_args={"n_procs": 4}, bypass_check=True)

    # Check output vs ref
    subjects = [
        "sub-ADNI011S4105",
        "sub-ADNI023S4020",
        "sub-ADNI035S4082",
        "sub-ADNI128S4832",
    ]
    out_data_forward_def = [
        fspath(
            output_dir
            / "caps"
            / "subjects"
            / sub
            / "ses-M000"
            / "t1"
            / "spm"
            / "dartel"
            / "group-UnitTest"
            / (
                sub
                + "_ses-M000_T1w_target-UnitTest_transformation-forward_deformation.nii.gz"
            )
        )
        for sub in subjects
    ]
    ref_data_forward_def = [
        fspath(
            ref_dir
            / (
                sub
                + "_ses-M000_T1w_target-UnitTest_transformation-forward_deformation.nii.gz"
            )
        )
        for sub in subjects
    ]

    for i in range(len(out_data_forward_def)):
        assert likeliness_measure(
            out_data_forward_def[i], ref_data_forward_def[i], (1e-3, 0.25), (1e-2, 0.1)
        )


def run_T1VolumeParcellation(
    input_dir: Path, output_dir: Path, ref_dir: Path, working_dir: Path
) -> None:
    import shutil

    import numpy as np
    import pandas as pds

    from clinica.pipelines.t1_volume_parcellation.t1_volume_parcellation_pipeline import (
        T1VolumeParcellation,
    )

    # Copy necessary data from in to out
    shutil.copytree(input_dir / "caps", output_dir / "caps", copy_function=shutil.copy)

    # Instantiate pipeline
    parameters = {"group_label": "UnitTest"}
    pipeline = T1VolumeParcellation(
        caps_directory=fspath(output_dir / "caps"),
        tsv_file=fspath(input_dir / "subjects.tsv"),
        base_dir=fspath(working_dir),
        parameters=parameters,
    )
    pipeline.build()
    pipeline.run(plugin="MultiProc", plugin_args={"n_procs": 4}, bypass_check=True)

    out_files = [
        fspath(
            output_dir
            / "caps"
            / "subjects"
            / "sub-ADNI018S4696"
            / "ses-M000"
            / "t1"
            / "spm"
            / "dartel"
            / "group-UnitTest"
            / "atlas_statistics"
            / (
                "sub-ADNI018S4696_ses-M000_T1w_segm-graymatter_space-Ixi549Space_modulated-on_probability_space-"
                + atlas
                + "_map-graymatter_statistics.tsv"
            )
        )
        for atlas in pipeline.parameters["atlases"]
    ]
    ref_files = [
        fspath(
            ref_dir
            / (
                "sub-ADNI018S4696_ses-M000_T1w_segm-graymatter_space-Ixi549Space_modulated-on_probability_space-"
                + atlas
                + "_map-graymatter_statistics.tsv"
            )
        )
        for atlas in pipeline.parameters["atlases"]
    ]

    for i in range(len(out_files)):
        out_csv = pds.read_csv(out_files[i], sep="\t")
        ref_csv = pds.read_csv(ref_files[i], sep="\t")
        assert np.allclose(
            np.array(out_csv.mean_scalar),
            np.array(ref_csv.mean_scalar),
            rtol=1e-8,
            equal_nan=True,
        )


def run_T1Linear(
    input_dir: Path, output_dir: Path, ref_dir: Path, working_dir: Path
) -> None:
    from os.path import abspath, dirname, join

    from clinica.pipelines.t1_linear.anat_linear_pipeline import AnatLinear

    parameters = {"uncropped_image": False}
    # Instantiate pipeline
    pipeline = AnatLinear(
        bids_directory=fspath(input_dir / "bids"),
        caps_directory=fspath(output_dir / "caps"),
        tsv_file=fspath(input_dir / "subjects.tsv"),
        base_dir=fspath(working_dir),
        parameters=parameters,
        name="t1-linear",
    )
    pipeline.run(plugin="MultiProc", plugin_args={"n_procs": 4}, bypass_check=True)

    compare_folders(output_dir / "caps", ref_dir / "caps", output_dir)


def run_FlairLinear(
    input_dir: Path, output_dir: Path, ref_dir: Path, working_dir: Path
) -> None:
    from os.path import abspath, dirname, join

    from clinica.pipelines.t1_linear.anat_linear_pipeline import AnatLinear

    parameters = {"uncropped_image": False}
    # Instantiate pipeline
    pipeline = AnatLinear(
        bids_directory=fspath(input_dir / "bids"),
        caps_directory=fspath(output_dir / "caps"),
        base_dir=fspath(working_dir),
        parameters=parameters,
        name="flair-linear",
    )
    pipeline.run(plugin="MultiProc", plugin_args={"n_procs": 4}, bypass_check=True)

    compare_folders(output_dir / "caps", ref_dir / "caps", output_dir)


def run_T1FreeSurferTemplate(
    input_dir: Path, output_dir: Path, ref_dir: Path, working_dir: Path
) -> None:
    # Data for this functional test comes from https://openneuro.org/datasets/ds000204
    # sub-01 was duplicated into to sub-02 with one session in order to test the "one time point" case
    import shutil

    from clinica.pipelines.t1_freesurfer_longitudinal.t1_freesurfer_template_pipeline import (
        T1FreeSurferTemplate,
    )

    # Copy necessary data from in to out
    shutil.copytree(input_dir / "caps", output_dir / "caps", copy_function=shutil.copy)

    pipeline = T1FreeSurferTemplate(
        caps_directory=fspath(output_dir / "caps"),
        tsv_file=fspath(input_dir / "subjects.tsv"),
        base_dir=fspath(working_dir),
    )
    pipeline.run(plugin="MultiProc", plugin_args={"n_procs": 2}, bypass_check=True)

    # We only check that folders are the same meaning that FreeSurfer finished without error
    # surf/ folder is ignored because it contains sym links that makes hard to check with ref data
    # (sym links of ref data are ignored after rsync on CI machines)
    def path_to_caps_fs(part_id: str, long_id: str) -> Path:

        output_folder = (
            Path("caps")
            / "subjects"
            / part_id
            / long_id
            / "freesurfer_unbiased_template"
        )
        return output_folder

    for (p_id, l_id) in zip(["sub-01", "sub-02"], ["long-20112015", "long-2011"]):

        folder1 = path_to_caps_fs(p_id, l_id) / (p_id + "_" + l_id)
        compare_folders(
            output_dir / folder1 / "label",
            ref_dir / folder1 / "label",
            output_dir,
        )
        compare_folders(
            output_dir / folder1 / "mri",
            ref_dir / folder1 / "mri",
            output_dir,
        )
        compare_folders(
            output_dir / folder1 / "stats",
            ref_dir / folder1 / "stats",
            output_dir,
        )


def run_T1FreeSurferLongitudinalCorrection(
    input_dir: Path, output_dir: Path, ref_dir: Path, working_dir: Path
) -> None:
    # Data for this functional test comes from https://openneuro.org/datasets/ds000204
    import shutil

    from clinica.pipelines.t1_freesurfer_longitudinal.t1_freesurfer_longitudinal_correction_pipeline import (
        T1FreeSurferLongitudinalCorrection,
    )

    # Copy necessary data from in to out
    shutil.copytree(input_dir / "caps", output_dir / "caps", copy_function=shutil.copy)

    pipeline = T1FreeSurferLongitudinalCorrection(
        caps_directory=fspath(output_dir / "caps"),
        tsv_file=fspath(input_dir / "subjects.tsv"),
        base_dir=fspath(working_dir),
    )
    pipeline.run(bypass_check=True)

    # We only check that folders are the same meaning that FreeSurfer finished without error
    # surf/ folder is ignored because it contains sym links that makes hard to check with ref data
    # (sym links of ref data are ignored after rsync on CI machines)
    def path_to_caps_fs(part_id: str, sess_id: str, long_id: str) -> Path:

        output_folder = (
            Path("caps")
            / "subjects"
            / part_id
            / sess_id
            / "t1"
            / long_id
            / "freesurfer_longitudinal"
        )
        return output_folder

    folder1 = path_to_caps_fs("sub-01", "ses-2011", "long-20112015")
    compare_folders(
        output_dir / folder1 / "regional_measures",
        ref_dir / folder1 / "regional_measures",
        output_dir,
    )
    compare_folders(
        output_dir / folder1 / "sub-01_ses-2011.long.sub-01_long-20112015" / "label",
        ref_dir / folder1 / "sub-01_ses-2011.long.sub-01_long-20112015" / "label",
        output_dir,
    )
    compare_folders(
        output_dir / folder1 / "sub-01_ses-2011.long.sub-01_long-20112015" / "mri",
        ref_dir / folder1 / "sub-01_ses-2011.long.sub-01_long-20112015" / "mri",
        output_dir,
    )
    compare_folders(
        output_dir / folder1 / "sub-01_ses-2011.long.sub-01_long-20112015" / "stats",
        ref_dir / folder1 / "sub-01_ses-2011.long.sub-01_long-20112015" / "stats",
        output_dir,
    )
