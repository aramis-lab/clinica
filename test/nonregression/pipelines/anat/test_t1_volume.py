import shutil
from os import fspath
from pathlib import Path
from test.nonregression.testing_tools import configure_paths, likeliness_measure
from typing import Tuple

import numpy as np
import pandas as pd


def test_t1_volume_tissue_segmentation(cmdopt, tmp_path):
    base_dir = Path(cmdopt["input"])
    working_dir = Path(cmdopt["wd"])
    input_dir, tmp_dir, ref_dir = configure_paths(
        base_dir, tmp_path, "T1VolumeTissueSegmentation"
    )
    run_t1_volume_tissue_segmentation(input_dir, tmp_dir, ref_dir, working_dir)


def run_t1_volume_tissue_segmentation(
    input_dir: Path, output_dir: Path, ref_dir: Path, working_dir: Path
) -> None:
    from clinica.pipelines.t1_volume_tissue_segmentation.t1_volume_tissue_segmentation_pipeline import (
        T1VolumeTissueSegmentation,
    )

    pipeline = T1VolumeTissueSegmentation(
        bids_directory=fspath(input_dir / "bids"),
        caps_directory=fspath(output_dir / "caps"),
        tsv_file=fspath(input_dir / "subjects.tsv"),
        parameters={"skip_question": False},
        base_dir=fspath(working_dir),
    )
    pipeline.build()
    pipeline.run(bypass_check=True)

    out_file = fspath(
        output_dir
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


def test_t1_volume_create_dartel(cmdopt, tmp_path):
    base_dir = Path(cmdopt["input"])
    working_dir = Path(cmdopt["wd"])
    input_dir, tmp_dir, ref_dir = configure_paths(
        base_dir, tmp_path, "T1VolumeCreateDartel"
    )
    run_t1_volume_create_dartel(input_dir, tmp_dir, ref_dir, working_dir)


def run_t1_volume_create_dartel(
    input_dir: Path, output_dir: Path, ref_dir: Path, working_dir: Path
) -> None:
    from clinica.pipelines.t1_volume_create_dartel.t1_volume_create_dartel_pipeline import (
        T1VolumeCreateDartel,
    )

    shutil.copytree(input_dir / "caps", output_dir / "caps", copy_function=shutil.copy)

    pipeline = T1VolumeCreateDartel(
        bids_directory=fspath(input_dir / "bids"),
        caps_directory=fspath(output_dir / "caps"),
        tsv_file=fspath(input_dir / "subjects.tsv"),
        base_dir=fspath(working_dir),
        group_label="UnitTest",
    )
    pipeline.build()
    pipeline.run(plugin="MultiProc", plugin_args={"n_procs": 4}, bypass_check=True)

    out_template = fspath(
        output_dir / "caps/groups/group-UnitTest/t1/group-UnitTest_template.nii.gz"
    )
    ref_template = fspath(ref_dir / "group-UnitTest_template.nii.gz")
    assert likeliness_measure(out_template, ref_template, (1e-3, 0.1), (1e-2, 0.1))

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
        for sub in _get_subjects()
    ]
    ref_data_forward_def = [
        fspath(
            ref_dir
            / (
                sub
                + "_ses-M000_T1w_target-UnitTest_transformation-forward_deformation.nii.gz"
            )
        )
        for sub in _get_subjects()
    ]

    for out, ref in zip(out_data_forward_def, ref_data_forward_def):
        assert likeliness_measure(out, ref, (1e-3, 0.25), (1e-2, 0.1))


def _get_subjects() -> Tuple[str, ...]:
    return (
        "sub-ADNI011S4105",
        "sub-ADNI023S4020",
        "sub-ADNI035S4082",
        "sub-ADNI128S4832",
    )


def test_t1_volume_dartel_to_mni(cmdopt, tmp_path):
    base_dir = Path(cmdopt["input"])
    working_dir = Path(cmdopt["wd"])
    input_dir, tmp_dir, ref_dir = configure_paths(
        base_dir, tmp_path, "T1VolumeDartel2MNI"
    )
    run_t1_volume_dartel_to_mni(input_dir, tmp_dir, ref_dir, working_dir)


def run_t1_volume_dartel_to_mni(
    input_dir: Path, output_dir: Path, ref_dir: Path, working_dir: Path
) -> None:
    from clinica.pipelines.t1_volume_dartel2mni.t1_volume_dartel2mni_pipeline import (
        T1VolumeDartel2MNI,
    )

    shutil.copytree(input_dir / "caps", output_dir / "caps", copy_function=shutil.copy)

    pipeline = T1VolumeDartel2MNI(
        bids_directory=fspath(input_dir / "bids"),
        caps_directory=fspath(output_dir / "caps"),
        tsv_file=fspath(input_dir / "subjects.tsv"),
        base_dir=fspath(working_dir),
        group_label="UnitTest",
    )
    pipeline.build()
    pipeline.run(plugin="MultiProc", plugin_args={"n_procs": 4}, bypass_check=True)

    out_data_gm_mni = [
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
        for sub in _get_subjects()
    ]
    ref_data_gm_mni = [
        fspath(
            ref_dir
            / (
                sub
                + "_ses-M000_T1w_segm-graymatter_space-Ixi549Space_modulated-on_fwhm-8mm_probability.nii.gz"
            )
        )
        for sub in _get_subjects()
    ]
    for out, ref in zip(out_data_gm_mni, ref_data_gm_mni):
        assert likeliness_measure(out, ref, (1e-4, 0.15), (1, 0.02))


def test_t1_volume_register_dartel(cmdopt, tmp_path):
    base_dir = Path(cmdopt["input"])
    working_dir = Path(cmdopt["wd"])
    input_dir, tmp_dir, ref_dir = configure_paths(
        base_dir, tmp_path, "T1VolumeRegisterDartel"
    )
    run_t1_volume_register_dartel(input_dir, tmp_dir, ref_dir, working_dir)


def run_t1_volume_register_dartel(
    input_dir: Path, output_dir: Path, ref_dir: Path, working_dir: Path
) -> None:
    from clinica.pipelines.t1_volume_register_dartel.t1_volume_register_dartel_pipeline import (
        T1VolumeRegisterDartel,
    )

    shutil.copytree(input_dir / "caps", output_dir / "caps", copy_function=shutil.copy)

    pipeline = T1VolumeRegisterDartel(
        bids_directory=fspath(input_dir / "bids"),
        caps_directory=fspath(output_dir / "caps"),
        tsv_file=fspath(input_dir / "subjects.tsv"),
        base_dir=fspath(working_dir),
        group_label="UnitTest",
    )
    pipeline.build()
    pipeline.run(plugin="MultiProc", plugin_args={"n_procs": 4}, bypass_check=True)

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
        for sub in _get_subjects()
    ]
    ref_data_forward_def = [
        fspath(
            ref_dir
            / (
                sub
                + "_ses-M000_T1w_target-UnitTest_transformation-forward_deformation.nii.gz"
            )
        )
        for sub in _get_subjects()
    ]

    for out, ref in zip(out_data_forward_def, ref_data_forward_def):
        assert likeliness_measure(out, ref, (1e-3, 0.25), (1e-2, 0.1))


def test_t1_volume_parcellation(cmdopt, tmp_path):
    base_dir = Path(cmdopt["input"])
    working_dir = Path(cmdopt["wd"])
    input_dir, tmp_dir, ref_dir = configure_paths(
        base_dir, tmp_path, "T1VolumeParcellation"
    )
    run_t1_volume_parcellation(input_dir, tmp_dir, ref_dir, working_dir)


def run_t1_volume_parcellation(
    input_dir: Path, output_dir: Path, ref_dir: Path, working_dir: Path
) -> None:
    from clinica.pipelines.t1_volume_parcellation.t1_volume_parcellation_pipeline import (
        T1VolumeParcellation,
    )

    shutil.copytree(input_dir / "caps", output_dir / "caps", copy_function=shutil.copy)

    pipeline = T1VolumeParcellation(
        caps_directory=fspath(output_dir / "caps"),
        tsv_file=fspath(input_dir / "subjects.tsv"),
        base_dir=fspath(working_dir),
        group_label="UnitTest",
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

    for out, ref in zip(out_files, ref_files):
        out_csv = pd.read_csv(out, sep="\t")
        ref_csv = pd.read_csv(ref, sep="\t")
        assert np.allclose(
            np.array(out_csv.mean_scalar),
            np.array(ref_csv.mean_scalar),
            rtol=1e-8,
            equal_nan=True,
        )
