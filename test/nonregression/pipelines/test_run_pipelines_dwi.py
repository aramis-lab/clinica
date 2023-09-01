# coding: utf8

"""
This file contains a set of functional tests designed to check the correct execution of the pipeline and the
different functions available in Clinica
"""

import warnings
from os import fspath
from pathlib import Path
from test.nonregression.testing_tools import configure_paths, similarity_measure

import nibabel as nib
import numpy as np
import pytest
from numpy.testing import assert_array_almost_equal

from clinica.utils.testing_utils import (
    assert_large_nifti_almost_equal,
    assert_nifti_almost_equal,
)

# Determine location for working_directory
warnings.filterwarnings("ignore")


@pytest.mark.fast
def test_dwi_b0_flirt(cmdopt, tmp_path):
    from clinica.pipelines.dwi_preprocessing_using_t1.dwi_preprocessing_using_t1_workflows import (
        b0_flirt_pipeline,
    )

    base_dir = Path(cmdopt["input"])
    input_dir, tmp_dir, ref_dir = configure_paths(base_dir, tmp_path, "DWIB0Flirt")
    b0_flirt = b0_flirt_pipeline(num_b0s=11)
    b0_flirt.inputs.inputnode.in_file = str(input_dir / "sub-01_ses-M000_dwi_b0.nii.gz")
    (tmp_path / "tmp").mkdir()
    b0_flirt.base_dir = str(tmp_path / "tmp")
    b0_flirt.run()

    out_file = fspath(
        tmp_path
        / "tmp"
        / "b0_coregistration"
        / "concat_ref_moving"
        / "merged_files.nii.gz"
    )
    ref_file = fspath(ref_dir / "merged_files.nii.gz")

    assert similarity_measure(out_file, ref_file, 0.99)


@pytest.mark.slow
def test_dwi_epi_pipeline(cmdopt, tmp_path):
    from clinica.pipelines.dwi_preprocessing_using_t1.dwi_preprocessing_using_t1_workflows import (
        epi_pipeline,
    )

    base_dir = Path(cmdopt["input"])
    input_dir, tmp_dir, ref_dir = configure_paths(base_dir, tmp_path, "DWIEPIPipeline")
    (tmp_path / "tmp").mkdir()
    epi = epi_pipeline(
        base_dir=str(base_dir),
        output_dir=str(tmp_path / "tmp"),
        ants_random_seed=42,
        use_double_precision=False,
        delete_cache=True,
    )
    epi.inputs.inputnode.t1_filename = str(input_dir / "sub-01_ses-M000_T1w.nii.gz")
    epi.inputs.inputnode.dwi_filename = str(input_dir / "sub-01_ses-M000_dwi.nii.gz")
    epi.inputs.inputnode.b_vectors_filename = str(
        input_dir / "sub-01_ses-M000_dwi.bvec"
    )

    epi.run()

    out_file = fspath(
        tmp_path / "tmp" / "rotated_b_vectors" / "sub-01_ses-M000_dwi_rotated.bvec"
    )
    ref_file = fspath(ref_dir / "sub-01_ses-M000_dwi_rotated.bvec")
    out_bvecs = np.loadtxt(out_file)
    ref_bvecs = np.loadtxt(ref_file)

    assert_array_almost_equal(out_bvecs, ref_bvecs, decimal=2)

    out_file = fspath(
        tmp_path
        / "tmp"
        / "epi_corrected_dwi_image"
        / "Jacobian_image_maths_thresh_merged.nii.gz"
    )
    ref_file = fspath(ref_dir / "Jacobian_image_maths_thresh_merged.nii.gz")

    assert_large_nifti_almost_equal(out_file, ref_file)


@pytest.mark.slow
def test_dwi_perform_ants_registration(cmdopt, tmp_path):
    from test.nonregression.testing_tools import similarity_measure

    import nibabel as nib

    from clinica.pipelines.dwi_preprocessing_using_t1.dwi_preprocessing_using_t1_workflows import (
        perform_ants_registration,
    )

    base_dir = Path(cmdopt["input"])
    input_dir, tmp_dir, ref_dir = configure_paths(
        base_dir, tmp_path, "DWIANTSRegistration"
    )
    (tmp_path / "tmp").mkdir()
    ants_registration = perform_ants_registration(
        base_dir=str(base_dir),
        output_dir=str(tmp_path / "tmp"),
        ants_random_seed=42,  # Set the random seed to avoid stochastic results (requires ants >= 2.3.0)
    )
    ants_registration.inputs.inputnode.t1_filename = str(
        input_dir / "sub-01_ses-M000_T1w.nii.gz"
    )
    ants_registration.inputs.inputnode.dwi_filename = str(
        input_dir / "sub-01_ses-M000_dwi.nii.gz"
    )
    ants_registration.inputs.inputnode.b_vectors_filename = str(
        input_dir / "sub-01_ses-M000_dwi.bvec"
    )

    ants_registration.run()

    out_file = fspath(
        tmp_path / "tmp" / "epi_correction_deformation_field" / "transform1Warp.nii.gz"
    )
    ref_file = fspath(ref_dir / "transform1Warp.nii.gz")

    assert nib.load(out_file).shape == nib.load(ref_file).shape

    out_file = fspath(
        tmp_path / "tmp" / "epi_correction_image_warped" / "transformWarped.nii.gz"
    )
    ref_file = fspath(ref_dir / "transformWarped.nii.gz")

    assert similarity_measure(out_file, ref_file, 0.9)

    out_file = fspath(tmp_path / "tmp" / "merged_transforms" / "transform1Warp.nii.gz")
    ref_file = fspath(ref_dir / "merged_transform.nii.gz")

    assert nib.load(out_file).shape == nib.load(ref_file).shape

    out_file = fspath(
        tmp_path / "tmp" / "rotated_b_vectors" / "sub-01_ses-M000_dwi_rotated.bvec"
    )
    ref_file = fspath(ref_dir / "sub-01_ses-M000_dwi_rotated.bvec")
    out_bvecs = np.loadtxt(out_file)
    ref_bvecs = np.loadtxt(ref_file)

    assert_array_almost_equal(out_bvecs, ref_bvecs, decimal=2)


@pytest.mark.slow
def test_dwi_perform_dwi_epi_correction(cmdopt, tmp_path):
    from clinica.pipelines.dwi_preprocessing_using_t1.dwi_preprocessing_using_t1_workflows import (
        perform_dwi_epi_correction,
    )

    base_dir = Path(cmdopt["input"])
    input_dir, tmp_dir, ref_dir = configure_paths(
        base_dir, tmp_path, "DWIEPICorrection"
    )
    (tmp_path / "tmp").mkdir()
    epi_correction = perform_dwi_epi_correction(
        base_dir=str(base_dir),
        output_dir=str(tmp_path / "tmp"),
        delete_cache=True,
        use_double_precision=False,
    )
    epi_correction.inputs.inputnode.t1_filename = str(
        input_dir / "sub-01_ses-M000_T1w.nii.gz"
    )
    epi_correction.inputs.inputnode.dwi_filename = str(
        input_dir / "sub-01_ses-M000_dwi.nii.gz"
    )
    epi_correction.inputs.inputnode.merged_transforms = str(
        input_dir / "merged_transforms.nii.gz"
    )

    epi_correction.run(plugin="MultiProc", plugin_args={"n_procs": 8})

    out_file = fspath(
        tmp_path
        / "tmp"
        / "epi_corrected_dwi_image"
        / "Jacobian_image_maths_thresh_merged.nii.gz"
    )
    ref_file = fspath(ref_dir / "Jacobian_image_maths_thresh_merged.nii.gz")

    assert_large_nifti_almost_equal(out_file, ref_file)


@pytest.mark.slow
def test_dwi_eddy_fsl(cmdopt, tmp_path):
    from clinica.pipelines.dwi_preprocessing_using_t1.dwi_preprocessing_using_t1_workflows import (
        eddy_fsl_pipeline,
    )

    base_dir = Path(cmdopt["input"])
    input_dir, tmp_dir, ref_dir = configure_paths(base_dir, tmp_path, "DWIEddyFSL")
    (tmp_path / "tmp").mkdir()
    eddy_fsl = eddy_fsl_pipeline(
        use_cuda=False,
        initrand=True,
        compute_mask=True,
        output_dir=str(tmp_path / "tmp"),
    )
    eddy_fsl.inputs.inputnode.total_readout_time = 0.0342002
    eddy_fsl.inputs.inputnode.phase_encoding_direction = "y-"
    eddy_fsl.inputs.inputnode.dwi_filename = str(
        input_dir / "sub-01_ses-M000_dwi.nii.gz"
    )
    eddy_fsl.inputs.inputnode.b_values_filename = str(
        input_dir / "sub-01_ses-M000_dwi.bval"
    )
    eddy_fsl.inputs.inputnode.b_vectors_filename = str(
        input_dir / "sub-01_ses-M000_dwi.bvec"
    )
    eddy_fsl.inputs.inputnode.reference_b0 = str(
        input_dir / "sub-01_ses-M000_dwi_b0.nii.gz"
    )

    eddy_fsl.run()

    out_file = fspath(tmp_path / "tmp" / "out_corrected" / "eddy_corrected.nii.gz")
    ref_file = fspath(ref_dir / "eddy_corrected.nii.gz")

    assert similarity_measure(out_file, ref_file, 0.97)

    out_file = fspath(
        tmp_path / "tmp" / "out_rotated_bvecs" / "eddy_corrected.eddy_rotated_bvecs"
    )
    ref_file = fspath(ref_dir / "eddy_corrected.eddy_rotated_bvecs")
    out_bvecs = np.loadtxt(out_file)
    ref_bvecs = np.loadtxt(ref_file)

    assert_array_almost_equal(out_bvecs, ref_bvecs)


@pytest.mark.slow
def test_dwi_preprocessing_using_t1(cmdopt, tmp_path):
    base_dir = Path(cmdopt["input"])
    working_dir = Path(cmdopt["wd"])
    input_dir, tmp_dir, ref_dir = configure_paths(
        base_dir, tmp_path, "DWIPreprocessingUsingT1"
    )
    run_dwi_preprocessing_using_t1(input_dir, tmp_dir, ref_dir, working_dir)


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
    from clinica.pipelines.dwi_preprocessing_using_fmap.dwi_preprocessing_using_phasediff_fmap_workflows import (
        compute_reference_b0,
    )
    from clinica.utils.dwi import bids_dir_to_fsl_dir
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
    phase_encoding_direction = bids_dir_to_fsl_dir(phase_encoding_direction)

    wf = compute_reference_b0(
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
        out_file = fspath(tmp_path / "tmp" / folder / filename)
        ref_file = fspath(ref_dir / folder / filename)

        assert similarity_measure(out_file, ref_file, 0.99)


@pytest.mark.fast
def test_prepare_phasediff_fmap(cmdopt, tmp_path):
    """Test the pipeline prepare_phasediff_fmap which is part of
     step 2 of the pipeline DWIPreprocessingUsingPhaseDiff.

    This is a fast test which should run in less than a minute.
    """
    from clinica.pipelines.dwi_preprocessing_using_fmap.dwi_preprocessing_using_phasediff_fmap_workflows import (
        prepare_phasediff_fmap,
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

    wf = prepare_phasediff_fmap(output_dir=str(tmp_path / "tmp"))
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
    out_file = fspath(tmp_path / "tmp" / "calibrated_fmap" / out_filename)
    ref_file = fspath(ref_dir / "calibrated_fmap" / out_filename)

    assert similarity_measure(out_file, ref_file, 0.99)


@pytest.mark.fast
def test_dwi_calibrate_and_register_fmap(cmdopt, tmp_path):
    """Test step 2 of pipeline DWIPreprocessingUsingPhaseDiff.

    This is a fast test which should run in about 1 minute.
    """
    from clinica.pipelines.dwi_preprocessing_using_fmap.dwi_preprocessing_using_phasediff_fmap_workflows import (
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

    wf = calibrate_and_register_fmap(output_dir=str(tmp_path / "tmp"))
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
        out_file = fspath(tmp_path / "tmp" / folder / filename)
        ref_file = fspath(ref_dir / folder / filename)

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


def run_dwi_preprocessing_using_t1(
    input_dir: Path, output_dir: Path, ref_dir: Path, working_dir: Path
) -> None:
    from clinica.pipelines.dwi_preprocessing_using_t1.dwi_preprocessing_using_t1_pipeline import (
        DwiPreprocessingUsingT1,
    )

    caps_dir = output_dir / "caps"
    tsv = input_dir / "subjects.tsv"

    parameters = {
        "initrand": True,
        "low_bval": 5,
        "use_cuda": False,
    }
    pipeline = DwiPreprocessingUsingT1(
        bids_directory=fspath(input_dir / "bids"),
        caps_directory=fspath(caps_dir),
        tsv_file=fspath(tsv),
        base_dir=fspath(working_dir),
        parameters=parameters,
    )
    pipeline.build()
    pipeline.run(plugin="MultiProc", plugin_args={"n_procs": 4}, bypass_check=True)

    # Assert
    out_file = fspath(
        caps_dir
        / "subjects"
        / "sub-PREVDEMALS0010025PG"
        / "ses-M000"
        / "dwi"
        / "preprocessing"
        / "sub-PREVDEMALS0010025PG-ses-M000_dwi_space-T1w_preproc.nii.gz"
    )
    ref_file = fspath(
        ref_dir / "sub-PREVDEMALS0010025PG_ses-M000_dwi_space-T1w_preproc.nii.gz"
    )

    assert similarity_measure(out_file, ref_file, 0.97)


def run_dwi_preprocessing_using_phase_diff_field_map(
    input_dir: Path, output_dir: Path, ref_dir: Path, working_dir: Path
) -> None:
    import warnings

    from clinica.pipelines.dwi_preprocessing_using_fmap.dwi_preprocessing_using_phasediff_fmap_pipeline import (
        DwiPreprocessingUsingPhaseDiffFMap,
    )

    warnings.filterwarnings("ignore")

    caps_dir = output_dir / "caps"
    tsv = input_dir / "subjects.tsv"

    parameters = {
        "initrand": True,
        "low_bval": 5,
        "use_cuda": False,
        "delete_cache": True,
    }
    pipeline = DwiPreprocessingUsingPhaseDiffFMap(
        bids_directory=fspath(input_dir / "bids"),
        caps_directory=fspath(caps_dir),
        tsv_file=fspath(tsv),
        base_dir=fspath(working_dir),
        parameters=parameters,
    )
    pipeline.build()
    pipeline.run(plugin="MultiProc", plugin_args={"n_procs": 4}, bypass_check=True)

    # Assert
    out_file = fspath(
        caps_dir
        / "subjects"
        / "sub-PREVDEMALS0010025PG"
        / "ses-M000"
        / "dwi"
        / "preprocessing"
        / "sub-PREVDEMALS0010025PG_ses-M000_dwi_space-b0_preproc.nii.gz"
    )
    ref_file = fspath(
        ref_dir / "sub-PREVDEMALS0010025PG_ses-M000_dwi_space-b0_preproc.nii.gz"
    )

    assert similarity_measure(out_file, ref_file, 0.95)


def run_dwi_dti(
    input_dir: Path, output_dir: Path, ref_dir: Path, working_dir: Path
) -> None:
    import shutil

    import numpy as np
    import pandas as pds

    from clinica.pipelines.dwi_dti.dwi_dti_pipeline import DwiDti

    caps_dir = output_dir / "caps"
    tsv = input_dir / "subjects.tsv"

    # Copy necessary data from in to out
    shutil.copytree(input_dir / "caps", caps_dir, copy_function=shutil.copy)

    pipeline = DwiDti(
        caps_directory=fspath(caps_dir),
        tsv_file=fspath(tsv),
        base_dir=fspath(working_dir),
    )
    pipeline.build()
    pipeline.run(plugin="MultiProc", plugin_args={"n_procs": 4}, bypass_check=True)

    # Check files
    subject_id = "sub-PREVDEMALS0010025PG"
    maps = ["AD", "FA", "MD", "RD"]
    out_files = [
        fspath(
            caps_dir
            / "subjects"
            / subject_id
            / "ses-M000"
            / "dwi"
            / "dti_based_processing"
            / "atlas_statistics"
            / (
                subject_id
                + "_ses-M000_dwi_space-JHUDTI81_res-1x1x1_map-"
                + m
                + "_statistics.tsv"
            )
        )
        for m in maps
    ]
    ref_files = [
        fspath(
            ref_dir
            / (
                subject_id
                + "_ses-M000_dwi_space-JHUDTI81_res-1x1x1_map-"
                + m
                + "_statistics.tsv"
            )
        )
        for m in maps
    ]

    for i in range(len(out_files)):
        out_csv = pds.read_csv(out_files[i], sep="\t")
        out_mean_scalar = np.array(out_csv.mean_scalar)
        ref_csv = pds.read_csv(ref_files[i], sep="\t")
        ref_mean_scalar = np.array(ref_csv.mean_scalar)

        assert np.allclose(out_mean_scalar, ref_mean_scalar, rtol=0.025, equal_nan=True)


def run_dwi_connectome(
    input_dir: Path, output_dir: Path, ref_dir: Path, working_dir: Path
) -> None:
    import shutil

    from clinica.pipelines.dwi_connectome.dwi_connectome_pipeline import DwiConnectome

    caps_dir = output_dir / "caps"
    tsv = input_dir / "subjects.tsv"

    # Copy necessary data from in to out
    shutil.copytree(input_dir / "caps", caps_dir, copy_function=shutil.copy)

    parameters = {"n_tracks": 1000}
    pipeline = DwiConnectome(
        caps_directory=fspath(caps_dir),
        tsv_file=fspath(tsv),
        base_dir=fspath(working_dir),
        parameters=parameters,
    )
    pipeline.build()
    pipeline.run(plugin="MultiProc", plugin_args={"n_procs": 4}, bypass_check=True)

    # Check files
    atlases = ["desikan", "destrieux"]
    subject_id = "sub-PREVDEMALS0010025PG"
    session_id = "ses-M000"

    out_fod_file = fspath(
        caps_dir
        / "subjects"
        / subject_id
        / session_id
        / "dwi"
        / "connectome_based_processing"
        / (subject_id + "_" + session_id + "_dwi_space-b0_model-CSD_diffmodel.nii.gz")
    )
    ref_fod_file = fspath(
        ref_dir
        / (subject_id + "_" + session_id + "_dwi_space-b0_model-CSD_diffmodel.nii.gz"),
    )

    out_parc_files = [
        fspath(
            caps_dir
            / "subjects"
            / subject_id
            / session_id
            / "dwi"
            / "connectome_based_processing"
            / (
                subject_id
                + "_"
                + session_id
                + "_dwi_space-b0_atlas-"
                + a
                + "_parcellation.nii.gz"
            )
        )
        for a in atlases
    ]
    ref_parc_files = [
        fspath(
            ref_dir
            / (
                subject_id
                + "_"
                + session_id
                + "_dwi_space-b0_atlas-"
                + a
                + "_parcellation.nii.gz"
            )
        )
        for a in atlases
    ]

    assert similarity_measure(out_fod_file, ref_fod_file, 0.97)

    for i in range(len(out_parc_files)):
        assert similarity_measure(out_parc_files[i], ref_parc_files[i], 0.955)
