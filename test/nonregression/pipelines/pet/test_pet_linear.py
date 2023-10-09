import shutil
from os import fspath
from pathlib import Path
from test.nonregression.testing_tools import compare_folders, configure_paths

from clinica.utils.pet import Tracer


def test_pet_linear(cmdopt, tmp_path):
    base_dir = Path(cmdopt["input"])
    working_dir = Path(cmdopt["wd"])
    input_dir, tmp_dir, ref_dir = configure_paths(base_dir, tmp_path, "PETLinear")
    run_pet_linear(input_dir, tmp_dir, ref_dir, working_dir)


def run_pet_linear(
    input_dir: Path, output_dir: Path, ref_dir: Path, working_dir: Path
) -> None:
    from clinica.pipelines.pet_linear.pet_linear_pipeline import PETLinear

    shutil.copytree(input_dir / "caps", output_dir / "caps", copy_function=shutil.copy)

    parameters = {"acq_label": Tracer.FDG, "suvr_reference_region": "pons"}

    pipeline = PETLinear(
        bids_directory=fspath(input_dir / "bids"),
        caps_directory=fspath(output_dir / "caps"),
        tsv_file=fspath(input_dir / "subjects.tsv"),
        base_dir=fspath(working_dir),
        parameters=parameters,
    )
    pipeline.run(plugin="MultiProc", plugin_args={"n_procs": 4}, bypass_check=True)

    compare_folders(output_dir / "caps", ref_dir / "caps", output_dir)
