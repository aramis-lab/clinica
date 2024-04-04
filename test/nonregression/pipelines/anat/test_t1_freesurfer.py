import shutil
from os import fspath
from pathlib import Path
from test.nonregression.testing_tools import compare_folders, configure_paths

import pytest


@pytest.mark.slow
def test_t1_freesurfer_cross_sectional(cmdopt, tmp_path):
    base_dir = Path(cmdopt["input"])
    working_dir = Path(cmdopt["wd"])
    input_dir, tmp_dir, ref_dir = configure_paths(base_dir, tmp_path, "T1FreeSurfer")
    run_t1_freesurfer_cross_sectional(input_dir, tmp_dir, ref_dir, working_dir)


def run_t1_freesurfer_cross_sectional(
    input_dir: Path, output_dir: Path, ref_dir: Path, working_dir: Path
) -> None:
    """Run the T1FreeSurfer pipeline on test data.

    Notes
    -----
    Data for this functional test comes from https://openneuro.org/datasets/ds000204

    We only check that folders are the same meaning that FreeSurfer finished without error.
    surf/ folder is ignored because it contains symlinks that makes hard to check with ref data
    (symlinks of ref data are ignored after rsync on CI machines).
    """
    from clinica.pipelines.anatomical.freesurfer.t1.pipeline import T1FreeSurfer

    parameters = {"recon_all_args": "-qcache", "skip_question": False}

    pipeline = T1FreeSurfer(
        bids_directory=fspath(input_dir / "bids"),
        caps_directory=fspath(output_dir / "caps"),
        tsv_file=fspath(input_dir / "subjects.tsv"),
        parameters=parameters,
        base_dir=fspath(working_dir),
    )
    pipeline.run(bypass_check=True)

    folder = _get_path_to_caps_freesurfer_cross_sectional("sub-01", "ses-2011")
    compare_folders(
        output_dir / folder / "regional_measures",
        ref_dir / folder / "regional_measures",
        output_dir,
    )
    for sub_folder in ("label", "mri", "stats"):
        compare_folders(
            output_dir / folder / "sub-01_ses-2011" / sub_folder,
            ref_dir / folder / "sub-01_ses-2011" / sub_folder,
            output_dir,
        )


def _get_path_to_caps_freesurfer_cross_sectional(part_id: str, session_id: str) -> Path:
    return (
        Path("caps")
        / "subjects"
        / part_id
        / session_id
        / "t1"
        / "freesurfer_cross_sectional"
    )


@pytest.mark.slow
def test_t1_freesurfer_template(cmdopt, tmp_path):
    base_dir = Path(cmdopt["input"])
    working_dir = Path(cmdopt["wd"])
    input_dir, tmp_dir, ref_dir = configure_paths(
        base_dir, tmp_path, "T1FreeSurferTemplate"
    )
    run_t1_freesurfer_template(input_dir, tmp_dir, ref_dir, working_dir)


def run_t1_freesurfer_template(
    input_dir: Path, output_dir: Path, ref_dir: Path, working_dir: Path
) -> None:
    """Run the T1FreeSurferTemplate pipeline on test data.

    Notes
    -----
    Data for this functional test comes from https://openneuro.org/datasets/ds000204

    sub-01 was duplicated into to sub-02 with one session in order to test the "one time point" case.

    We only check that folders are the same meaning that FreeSurfer finished without error.
    surf/ folder is ignored because it contains symlinks that makes hard to check with ref data
    (symlinks of ref data are ignored after rsync on CI machines).
    """
    from clinica.pipelines.anatomical.freesurfer.longitudinal.template.pipeline import (
        T1FreeSurferTemplate,
    )

    shutil.copytree(input_dir / "caps", output_dir / "caps", copy_function=shutil.copy)

    pipeline = T1FreeSurferTemplate(
        caps_directory=fspath(output_dir / "caps"),
        tsv_file=fspath(input_dir / "subjects.tsv"),
        base_dir=fspath(working_dir),
    )
    pipeline.run(plugin="MultiProc", plugin_args={"n_procs": 2}, bypass_check=True)

    for part_id, long_id in zip(["sub-01", "sub-02"], ["long-20112015", "long-2011"]):
        folder = (
            _get_path_to_caps_freesurfer_template(part_id, long_id)
            / f"{part_id}_{long_id}"
        )
        for sub_folder in ("label", "mri", "stats"):
            compare_folders(
                output_dir / folder / sub_folder,
                ref_dir / folder / sub_folder,
                output_dir,
            )


def _get_path_to_caps_freesurfer_template(part_id: str, long_id: str) -> Path:
    return (
        Path("caps") / "subjects" / part_id / long_id / "freesurfer_unbiased_template"
    )


@pytest.mark.slow
def test_t1_freesurfer_longitudinal_correction(cmdopt, tmp_path):
    base_dir = Path(cmdopt["input"])
    working_dir = Path(cmdopt["wd"])
    input_dir, tmp_dir, ref_dir = configure_paths(
        base_dir, tmp_path, "T1FreeSurferLongitudinalCorrection"
    )
    run_t1_freesurfer_longitudinal_correction(input_dir, tmp_dir, ref_dir, working_dir)


def run_t1_freesurfer_longitudinal_correction(
    input_dir: Path, output_dir: Path, ref_dir: Path, working_dir: Path
) -> None:
    """Run the T1FreeSurferLongitudinalCorrection pipeline on test data.

    Notes
    -----
    Data for this functional test comes from https://openneuro.org/datasets/ds000204

    We only check that folders are the same meaning that FreeSurfer finished without error.
    surf/ folder is ignored because it contains symlinks that makes hard to check with ref data
    (symlinks of ref data are ignored after rsync on CI machines).
    """
    from clinica.pipelines.anatomical.freesurfer.longitudinal.correction.pipeline import (
        T1FreeSurferLongitudinalCorrection,
    )

    shutil.copytree(input_dir / "caps", output_dir / "caps", copy_function=shutil.copy)

    pipeline = T1FreeSurferLongitudinalCorrection(
        caps_directory=fspath(output_dir / "caps"),
        tsv_file=fspath(input_dir / "subjects.tsv"),
        base_dir=fspath(working_dir),
    )
    pipeline.run(bypass_check=True)

    folder = _get_path_to_caps_freesurfer_longitudinal(
        "sub-01", "ses-2011", "long-20112015"
    )
    compare_folders(
        output_dir / folder / "regional_measures",
        ref_dir / folder / "regional_measures",
        output_dir,
    )
    for sub_folder in ("label", "mri", "stats"):
        compare_folders(
            output_dir
            / folder
            / "sub-01_ses-2011.long.sub-01_long-20112015"
            / sub_folder,
            ref_dir / folder / "sub-01_ses-2011.long.sub-01_long-20112015" / sub_folder,
            output_dir,
        )


def _get_path_to_caps_freesurfer_longitudinal(
    part_id: str, session_id: str, long_id: str
) -> Path:
    return (
        Path("caps")
        / "subjects"
        / part_id
        / session_id
        / "t1"
        / long_id
        / "freesurfer_longitudinal"
    )
