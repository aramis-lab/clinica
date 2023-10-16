"""Non regression tests for sub-workflows of the DWIPreprocessingUsingT1 pipeline."""

from os import fspath
from pathlib import Path
from test.nonregression.testing_tools import (
    configure_paths,
    similarity_measure,
    similarity_measure_large_nifti,
)

import numpy as np
import pytest
from numpy.testing import assert_array_almost_equal


@pytest.mark.fast
def test_dwi_b0_flirt(cmdopt, tmp_path):
    from clinica.pipelines.dwi_preprocessing_using_t1.workflows import b0_flirt_pipeline

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
    from clinica.pipelines.dwi_preprocessing_using_t1.workflows import epi_pipeline

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

    out_b_vectors_file = fspath(
        tmp_path / "tmp" / "rotated_b_vectors" / "sub-01_ses-M000_dwi_rotated.bvec"
    )
    ref_b_vectors_file = fspath(ref_dir / "sub-01_ses-M000_dwi_rotated.bvec")
    out_b_vectors = np.loadtxt(out_b_vectors_file)
    ref_b_vectors = np.loadtxt(ref_b_vectors_file)

    assert_array_almost_equal(out_b_vectors, ref_b_vectors, decimal=2)

    out_jacobian = fspath(
        tmp_path
        / "tmp"
        / "epi_corrected_dwi_image"
        / "Jacobian_image_maths_thresh_merged.nii.gz"
    )
    ref_jacobian = fspath(ref_dir / "Jacobian_image_maths_thresh_merged.nii.gz")

    similarity_measure_large_nifti(out_jacobian, ref_jacobian, 0.95)


@pytest.mark.slow
def test_dwi_perform_ants_registration(cmdopt, tmp_path):
    from test.nonregression.testing_tools import similarity_measure

    import nibabel as nib

    from clinica.pipelines.dwi_preprocessing_using_t1.workflows import (
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

    out_deformation_field = fspath(
        tmp_path / "tmp" / "epi_correction_deformation_field" / "transform1Warp.nii.gz"
    )
    ref_deformation_field = fspath(ref_dir / "transform1Warp.nii.gz")

    assert (
        nib.load(out_deformation_field).shape == nib.load(ref_deformation_field).shape
    )

    out_warped_transform = fspath(
        tmp_path / "tmp" / "epi_correction_image_warped" / "transformWarped.nii.gz"
    )
    ref_warped_transform = fspath(ref_dir / "transformWarped.nii.gz")

    assert similarity_measure(out_warped_transform, ref_warped_transform, 0.9)

    out_merged_transforms = fspath(
        tmp_path / "tmp" / "merged_transforms" / "transform1Warp.nii.gz"
    )
    ref_merged_transforms = fspath(ref_dir / "merged_transform.nii.gz")

    assert (
        nib.load(out_merged_transforms).shape == nib.load(ref_merged_transforms).shape
    )

    out_b_vectors_file = fspath(
        tmp_path / "tmp" / "rotated_b_vectors" / "sub-01_ses-M000_dwi_rotated.bvec"
    )
    ref_b_vectors_file = fspath(ref_dir / "sub-01_ses-M000_dwi_rotated.bvec")
    out_b_vectors = np.loadtxt(out_b_vectors_file)
    ref_b_vectors = np.loadtxt(ref_b_vectors_file)

    assert_array_almost_equal(out_b_vectors, ref_b_vectors, decimal=2)


@pytest.mark.slow
def test_dwi_perform_dwi_epi_correction(cmdopt, tmp_path):
    from clinica.pipelines.dwi_preprocessing_using_t1.workflows import (
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

    similarity_measure_large_nifti(out_file, ref_file, 0.95)


@pytest.mark.slow
def test_dwi_eddy_fsl(cmdopt, tmp_path):
    from clinica.pipelines.dwi_preprocessing_using_t1.workflows import eddy_fsl_pipeline

    base_dir = Path(cmdopt["input"])
    input_dir, tmp_dir, ref_dir = configure_paths(base_dir, tmp_path, "DWIEddyFSL")
    (tmp_path / "tmp").mkdir()
    eddy_fsl = eddy_fsl_pipeline(
        base_dir=str(base_dir),
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
    out_b_vectors = np.loadtxt(out_file)
    ref_b_vectors = np.loadtxt(ref_file)

    assert_array_almost_equal(out_b_vectors, ref_b_vectors, decimal=3)
