from pathlib import Path
from test.nonregression.testing_tools import (
    configure_paths,
    similarity_measure,
)

import pytest


@pytest.mark.slow
def test_dwi_compute_reference_b0(cmdopt, tmp_path):
    """Test step 1 of pipeline DWIPreprocessingUsingPhaseDiff.

    This is a slow test. Expect approximately a 2.5 hours runtime.

    We use the same input files as for the DWIPreprocessingUsingPhaseDiff pipeline.

    However, since we are testing the step 1 in isolation, we don't have the
    BIDS querying logic available to us. The compute_reference_b0 pipeline is
    not BIDS-aware, it just expects three files:
        - the DWI image
        - the associated b-values
        - and the associated b-vectors

    These files are directly available in the input_dir (no BIDS structure, files
    are directly at the root level).

    It also needs the phase encoding direction and total readout time which are
    available in the JSON file (also available in the input_dir).
    In the main pipeline, the input node is performing the extraction of these
    metadata from the JSON file. Here, we perform this extraction directly with
    the `extract_metadata_from_json` and `bids_dir_to_fsl_dir`functions.

    The compute_reference_b0 workflow produces two outputs that are compared
    against their reference values:
        - The reference B0 volume
        - The brain mask computed on the B0 volume
    """
    from clinica.pipelines.dwi.preprocessing.fmap.workflows import (
        compute_reference_b0,
    )
    from clinica.pipelines.dwi.preprocessing.utils import _bids_dir_to_fsl_dir  # noqa
    from clinica.utils.filemanip import (
        extract_metadata_from_json,
        handle_missing_keys_dwi,
    )

    base_dir = Path(cmdopt["input"])
    input_dir, tmp_dir, ref_dir = configure_paths(
        base_dir, tmp_path, "DWIComputeReferenceB0"
    )
    (tmp_path / "tmp").mkdir()

    [total_readout_time, phase_encoding_direction] = extract_metadata_from_json(
        str(input_dir / "sub-01_ses-M000_dwi.json"),
        ["TotalReadoutTime", "PhaseEncodingDirection"],
        handle_missing_keys=handle_missing_keys_dwi,
    )
    phase_encoding_direction = _bids_dir_to_fsl_dir(phase_encoding_direction)

    wf = compute_reference_b0(
        base_dir=str(tmp_dir),
        b_value_threshold=5.0,
        use_cuda=False,
        initrand=False,
        output_dir=str(tmp_path / "tmp"),
        name="compute_reference_b0",
    )
    wf.inputs.inputnode.b_values_filename = str(input_dir / "sub-01_ses-M000_dwi.bval")
    wf.inputs.inputnode.b_vectors_filename = str(input_dir / "sub-01_ses-M000_dwi.bvec")
    wf.inputs.inputnode.dwi_filename = str(input_dir / "sub-01_ses-M000_dwi.nii.gz")
    wf.inputs.inputnode.image_id = "sub-01_ses-M000"
    wf.inputs.inputnode.total_readout_time = total_readout_time
    wf.inputs.inputnode.phase_encoding_direction = phase_encoding_direction

    wf.run()

    for folder, filename in zip(
        ["reference_b0", "brainmask"],
        ["sub-01_ses-M000_avg_b0_brain.nii.gz", "brainmask.nii.gz"],
    ):
        out_file = tmp_path / "tmp" / folder / filename
        ref_file = ref_dir / folder / filename

        assert similarity_measure(out_file, ref_file, 0.99)


@pytest.mark.fast
def test_prepare_phasediff_fmap(cmdopt, tmp_path):
    """Test the pipeline prepare_phasediff_fmap which is part of
     step 2 of the pipeline DWIPreprocessingUsingPhaseDiff.

    This is a fast test which should run in less than a minute.
    """
    from clinica.pipelines.dwi.preprocessing.fmap.workflows import (
        _prepare_phasediff_fmap,  # noqa
    )
    from clinica.utils.filemanip import extract_metadata_from_json

    base_dir = Path(cmdopt["input"])
    input_dir, tmp_dir, ref_dir = configure_paths(
        base_dir, tmp_path, "DWIPreparePhasediffFmap"
    )
    (tmp_path / "tmp").mkdir()

    [echo_time1, echo_time2] = extract_metadata_from_json(
        str(input_dir / "sub-01_ses-M000_phasediff.json"), ["EchoTime1", "EchoTime2"]
    )
    delta_echo_time = abs(echo_time2 - echo_time1)

    wf = _prepare_phasediff_fmap(
        base_dir=str(tmp_dir), output_dir=str(tmp_path / "tmp")
    )
    wf.inputs.input_node.fmap_mask = str(
        input_dir / "sub-01_ses-M000_magnitude1_corrected_brain_mask.nii.gz"
    )
    wf.inputs.input_node.fmap_phasediff = str(
        input_dir / "sub-01_ses-M000_phasediff.nii.gz"
    )
    wf.inputs.input_node.fmap_magnitude = str(
        input_dir / "sub-01_ses-M000_magnitude1_corrected_brain.nii.gz"
    )
    wf.inputs.input_node.delta_echo_time = str(delta_echo_time)

    wf.run()

    out_filename = (
        "sub-01_ses-M000_phasediff_rads_unwrapped_radsec_fieldmap_demean_maths.nii.gz"
    )
    out_file = tmp_path / "tmp" / "calibrated_fmap" / out_filename
    ref_file = ref_dir / "calibrated_fmap" / out_filename

    assert similarity_measure(out_file, ref_file, 0.99)


@pytest.mark.fast
def test_dwi_calibrate_and_register_fmap(cmdopt, tmp_path):
    """Test step 2 of pipeline DWIPreprocessingUsingPhaseDiff.

    This is a fast test which should run in about 1 minute.
    """
    from clinica.pipelines.dwi.preprocessing.fmap.workflows import (
        calibrate_and_register_fmap,
    )
    from clinica.utils.filemanip import extract_metadata_from_json

    base_dir = Path(cmdopt["input"])
    input_dir, tmp_dir, ref_dir = configure_paths(
        base_dir, tmp_path, "DWICalibrateAndRegisterFmap"
    )
    (tmp_path / "tmp").mkdir()

    [echo_time1, echo_time2] = extract_metadata_from_json(
        str(input_dir / "sub-01_ses-M000_phasediff.json"), ["EchoTime1", "EchoTime2"]
    )
    delta_echo_time = abs(echo_time2 - echo_time1)

    wf = calibrate_and_register_fmap(
        base_dir=str(tmp_dir), output_dir=str(tmp_path / "tmp")
    )
    wf.inputs.inputnode.reference_b0 = str(
        input_dir / "sub-01_ses-M000_avg_b0_brain.nii.gz"
    )
    wf.inputs.inputnode.bias_magnitude_fmap = str(
        input_dir / "sub-01_ses-M000_magnitude1.nii.gz"
    )
    wf.inputs.inputnode.fmap_phasediff = str(
        input_dir / "sub-01_ses-M000_phasediff.nii.gz"
    )
    wf.inputs.inputnode.delta_echo_time = delta_echo_time

    wf.run()

    for folder, filename in zip(
        [
            "smooth_calibrated_fmap",
            "bet_magnitude_fmap_registered_onto_b0",
            "registered_calibrated_fmap",
        ],
        [
            "sub-01_ses-M000_phasediff_rads_unwrapped_radsec_fieldmap_demean_maths_flirt_smooth.nii.gz",
            "sub-01_ses-M000_magnitude1_corrected_brain_flirt.nii.gz",
            "sub-01_ses-M000_phasediff_rads_unwrapped_radsec_fieldmap_demean_maths_flirt.nii.gz",
        ],
    ):
        out_file = tmp_path / "tmp" / folder / filename
        ref_file = ref_dir / folder / filename

        assert similarity_measure(out_file, ref_file, 0.98)


@pytest.mark.slow
def test_dwi_preprocessing_using_phase_diff_field_map(cmdopt, tmp_path):
    base_dir = Path(cmdopt["input"])
    working_dir = Path(cmdopt["wd"])
    input_dir, tmp_dir, ref_dir = configure_paths(
        base_dir,
        tmp_path,
        "DWIPreprocessingUsingPhaseDiffFieldmap",
    )
    run_dwi_preprocessing_using_phase_diff_field_map(
        input_dir, tmp_dir, ref_dir, working_dir
    )


def run_dwi_preprocessing_using_phase_diff_field_map(
    input_dir: Path, output_dir: Path, ref_dir: Path, working_dir: Path
) -> None:
    from clinica.pipelines.dwi.preprocessing.fmap.pipeline import (
        DwiPreprocessingUsingPhaseDiffFMap,
    )

    caps_dir = output_dir / "caps"
    parameters = {
        "initrand": True,
        "low_bval": 5,
        "use_cuda": False,
        "delete_cache": True,
    }
    pipeline = DwiPreprocessingUsingPhaseDiffFMap(
        bids_directory=str(input_dir / "bids"),
        caps_directory=str(caps_dir),
        tsv_file=str(input_dir / "subjects.tsv"),
        base_dir=str(working_dir),
        parameters=parameters,
    )
    pipeline.build()
    pipeline.run(plugin="MultiProc", plugin_args={"n_procs": 4}, bypass_check=True)

    out_file = (
        caps_dir
        / "subjects"
        / "sub-PREVDEMALS0010025PG"
        / "ses-M000"
        / "dwi"
        / "preprocessing"
        / "sub-PREVDEMALS0010025PG_ses-M000_dwi_space-b0_preproc.nii.gz"
    )
    ref_file = ref_dir / "sub-PREVDEMALS0010025PG_ses-M000_dwi_space-b0_preproc.nii.gz"

    assert similarity_measure(out_file, ref_file, 0.95)
