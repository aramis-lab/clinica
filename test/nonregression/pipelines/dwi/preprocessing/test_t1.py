from pathlib import Path
from test.nonregression.testing_tools import (
    configure_paths,
    similarity_measure,
    similarity_measure_large_nifti,
)

import nibabel as nib
import numpy as np
import pytest
from numpy.testing import assert_array_almost_equal


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

    out_file = (
        tmp_path
        / "tmp"
        / "b0_coregistration"
        / "concat_ref_moving"
        / "merged_files.nii.gz"
    )
    ref_file = ref_dir / "merged_files.nii.gz"

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

    out_file = (
        tmp_path / "tmp" / "rotated_b_vectors" / "sub-01_ses-M000_dwi_rotated.bvec"
    )
    ref_file = ref_dir / "sub-01_ses-M000_dwi_rotated.bvec"
    out_bvecs = np.loadtxt(out_file)
    ref_bvecs = np.loadtxt(ref_file)

    assert_array_almost_equal(out_bvecs, ref_bvecs, decimal=2)

    out_file = (
        tmp_path
        / "tmp"
        / "epi_corrected_dwi_image"
        / "Jacobian_image_maths_thresh_merged.nii.gz"
    )
    ref_file = ref_dir / "Jacobian_image_maths_thresh_merged.nii.gz"

    similarity_measure_large_nifti(out_file, ref_file, 0.95)


@pytest.mark.slow
def test_dwi_perform_ants_registration(cmdopt, tmp_path):
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

    out_file = (
        tmp_path / "tmp" / "epi_correction_deformation_field" / "transform1Warp.nii.gz"
    )
    ref_file = ref_dir / "transform1Warp.nii.gz"

    assert nib.load(out_file).shape == nib.load(ref_file).shape

    out_file = (
        tmp_path / "tmp" / "epi_correction_image_warped" / "transformWarped.nii.gz"
    )
    ref_file = ref_dir / "transformWarped.nii.gz"

    assert similarity_measure(out_file, ref_file, 0.9)

    out_file = tmp_path / "tmp" / "merged_transforms" / "transform1Warp.nii.gz"
    ref_file = ref_dir / "merged_transform.nii.gz"

    assert nib.load(out_file).shape == nib.load(ref_file).shape

    out_file = (
        tmp_path / "tmp" / "rotated_b_vectors" / "sub-01_ses-M000_dwi_rotated.bvec"
    )
    ref_file = ref_dir / "sub-01_ses-M000_dwi_rotated.bvec"
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

    out_file = (
        tmp_path
        / "tmp"
        / "epi_corrected_dwi_image"
        / "Jacobian_image_maths_thresh_merged.nii.gz"
    )
    ref_file = ref_dir / "Jacobian_image_maths_thresh_merged.nii.gz"

    similarity_measure_large_nifti(out_file, ref_file, 0.95)


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

    out_file = tmp_path / "tmp" / "out_corrected" / "eddy_corrected.nii.gz"
    ref_file = ref_dir / "eddy_corrected.nii.gz"

    assert similarity_measure(out_file, ref_file, 0.97)

    out_file = (
        tmp_path / "tmp" / "out_rotated_bvecs" / "eddy_corrected.eddy_rotated_bvecs"
    )
    ref_file = ref_dir / "eddy_corrected.eddy_rotated_bvecs"
    out_bvecs = np.loadtxt(out_file)
    ref_bvecs = np.loadtxt(ref_file)

    assert_array_almost_equal(out_bvecs, ref_bvecs, decimal=3)


@pytest.mark.slow
def test_dwi_preprocessing_using_t1(cmdopt, tmp_path):
    base_dir = Path(cmdopt["input"])
    working_dir = Path(cmdopt["wd"])
    input_dir, tmp_dir, ref_dir = configure_paths(
        base_dir, tmp_path, "DWIPreprocessingUsingT1"
    )
    run_dwi_preprocessing_using_t1(input_dir, tmp_dir, ref_dir, working_dir)


def run_dwi_preprocessing_using_t1(
    input_dir: Path, output_dir: Path, ref_dir: Path, working_dir: Path
) -> None:
    from clinica.pipelines.dwi_preprocessing_using_t1.dwi_preprocessing_using_t1_pipeline import (
        DwiPreprocessingUsingT1,
    )

    caps_dir = output_dir / "caps"
    parameters = {
        "initrand": True,
        "low_bval": 5,
        "use_cuda": False,
        "random_seed": 42,
    }
    pipeline = DwiPreprocessingUsingT1(
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
        / "sub-PREVDEMALS0010025PG_ses-M000_space-T1w_desc-preproc_dwi.nii.gz"
    )
    ref_file = (
        ref_dir / "sub-PREVDEMALS0010025PG_ses-M000_space-T1w_desc-preproc_dwi.nii.gz"
    )

    assert similarity_measure_large_nifti(out_file, ref_file, 0.99)
