import shutil
from pathlib import Path
from test.nonregression.testing_tools import (
    configure_paths,
    similarity_measure,
)

import numpy as np
import pandas as pd
import pytest
from numpy.testing import assert_array_almost_equal


@pytest.mark.slow
def test_dwi_dti(cmdopt, tmp_path):
    base_dir = Path(cmdopt["input"])
    working_dir = Path(cmdopt["wd"])
    input_dir, tmp_dir, ref_dir = configure_paths(base_dir, tmp_path, "DWIDTI")
    run_dwi_dti(input_dir, tmp_dir, ref_dir, working_dir)


@pytest.mark.slow
def test_dwi_connectome(cmdopt, tmp_path):
    base_dir = Path(cmdopt["input"])
    working_dir = Path(cmdopt["wd"])
    input_dir, tmp_dir, ref_dir = configure_paths(base_dir, tmp_path, "DWIConnectome")
    run_dwi_connectome(input_dir, tmp_dir, ref_dir, working_dir)


def run_dwi_dti(
    input_dir: Path, output_dir: Path, ref_dir: Path, working_dir: Path
) -> None:
    from clinica.pipelines.dwi_dti.pipeline import DwiDti
    from clinica.utils.bids import BIDSFileName
    from clinica.utils.dwi import DTIBasedMeasure

    caps_dir = output_dir / "caps"

    # Copy necessary data from in to out
    shutil.copytree(input_dir / "caps", caps_dir, copy_function=shutil.copy)

    pipeline = DwiDti(
        caps_directory=str(caps_dir),
        tsv_file=str(input_dir / "subjects.tsv"),
        base_dir=str(working_dir),
        parameters={"random_seed": 42},
    )
    pipeline.build()
    pipeline.run(plugin="MultiProc", plugin_args={"n_procs": 4}, bypass_check=True)

    filename = BIDSFileName.from_name(
        "sub-01_ses-M000_space-JHUDTI81_desc-preproc_res-1x1x1_statistics.tsv"
    )
    output = (
        caps_dir
        / "subjects"
        / "sub-01"
        / "ses-M000"
        / "dwi"
        / "dti_based_processing"
        / "atlas_statistics"
    )
    for measure in DTIBasedMeasure:
        filename.update_entity("map", measure.value)
        out_csv = pd.read_csv(output / filename.name, sep="\t")
        ref_csv = pd.read_csv(ref_dir / filename.name, sep="\t")
        assert_array_almost_equal(
            np.array(out_csv.mean_scalar), np.array(ref_csv.mean_scalar), decimal=2
        )


def run_dwi_connectome(
    input_dir: Path, output_dir: Path, ref_dir: Path, working_dir: Path
) -> None:
    from clinica.pipelines.dwi_connectome.pipeline import DwiConnectome
    from clinica.utils.bids import BIDSFileName

    caps_dir = output_dir / "caps"

    # Copy necessary data from in to out
    shutil.copytree(input_dir / "caps", caps_dir, copy_function=shutil.copy)

    parameters = {"n_tracks": 1000}
    pipeline = DwiConnectome(
        caps_directory=str(caps_dir),
        tsv_file=str(input_dir / "subjects.tsv"),
        base_dir=str(working_dir),
        parameters=parameters,
    )
    pipeline.build()
    pipeline.run(plugin="MultiProc", plugin_args={"n_procs": 4}, bypass_check=True)

    filename = BIDSFileName.from_name(
        "sub-01_ses-M000_space-b0_desc-preproc_model-CSD_diffmodel.nii.gz"
    )
    output_folder = (
        caps_dir
        / "subjects"
        / "sub-01"
        / "ses-M000"
        / "dwi"
        / "connectome_based_processing"
    )
    out_fod_file = output_folder / filename.name
    ref_fod_file = ref_dir / filename.name

    assert similarity_measure(out_fod_file, ref_fod_file, 0.97)

    for atlas in ("desikan", "destrieux"):
        filename.update_entity("atlas", atlas)
        filename.delete_entity("model")
        filename.suffix = "parcellation"
        assert similarity_measure(
            output_folder / filename.name,
            ref_dir / filename.name,
            0.955,
        )
