import shutil
from os import fspath
from pathlib import Path
from test.nonregression.testing_tools import configure_paths, likeliness_measure

from clinica.utils.pet import SUVRReferenceRegion, Tracer


def test_pet_volume(cmdopt, tmp_path):
    base_dir = Path(cmdopt["input"])
    working_dir = Path(cmdopt["wd"])
    input_dir, tmp_dir, ref_dir = configure_paths(base_dir, tmp_path, "PETVolume")
    run_pet_volume(input_dir, tmp_dir, ref_dir, working_dir)


def run_pet_volume(
    input_dir: Path, output_dir: Path, ref_dir: Path, working_dir: Path
) -> None:
    from clinica.pipelines.pet.volume.pipeline import PETVolume

    shutil.copytree(input_dir / "caps", output_dir / "caps", copy_function=shutil.copy)

    tracer = Tracer.FDG
    region = SUVRReferenceRegion.PONS

    parameters = {
        "group_label": "UnitTest",
        "acq_label": tracer,
        "suvr_reference_region": region,
        "skip_question": True,
    }
    pipeline = PETVolume(
        bids_directory=fspath(input_dir / "bids"),
        caps_directory=fspath(output_dir / "caps"),
        tsv_file=fspath(input_dir / "subjects.tsv"),
        base_dir=fspath(working_dir),
        parameters=parameters,
    )
    pipeline.build()
    pipeline.run(plugin="MultiProc", plugin_args={"n_procs": 4}, bypass_check=True)

    suffix = "_mask-brain_fwhm-8mm_pet.nii.gz"
    for subject in (
        "sub-ADNI011S4105",
        "sub-ADNI023S4020",
        "sub-ADNI035S4082",
        "sub-ADNI128S4832",
    ):
        output_folder = (
            output_dir
            / "caps"
            / "subjects"
            / subject
            / "ses-M000/pet/preprocessing/group-UnitTest"
        )
        assert likeliness_measure(
            output_folder
            / f"{subject}_ses-M000_trc-{tracer.value}_pet_space-Ixi549Space_suvr-{region.value}{suffix}",
            ref_dir
            / f"{subject}_ses-M000_trc-{tracer.value}_pet_space-Ixi549Space_suvr-{region.value}{suffix}",
            (1e-2, 0.25),
            (1e-1, 0.001),
        )
