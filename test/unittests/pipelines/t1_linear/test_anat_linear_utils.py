import pytest


@pytest.mark.parametrize("suffix", ["T1w", "FLAIR", "fooo"])
def test_get_substitutions_datasink(suffix):
    from clinica.pipelines.t1_linear.anat_linear_utils import (
        _get_substitutions_datasink,
    )

    bids_image_id = f"sub-ADNI022S0004_ses-M000_{suffix}"

    substitutions = _get_substitutions_datasink(bids_image_id, suffix)

    assert len(substitutions) == 3
    assert substitutions[0] == (
        f"sub-ADNI022S0004_ses-M000_{suffix}Warped_cropped.nii.gz",
        f"sub-ADNI022S0004_ses-M000_space-MNI152NLin2009cSym_desc-Crop_res-1x1x1_{suffix}.nii.gz",
    )
    assert substitutions[1] == (
        f"sub-ADNI022S0004_ses-M000_{suffix}0GenericAffine.mat",
        f"sub-ADNI022S0004_ses-M000_space-MNI152NLin2009cSym_res-1x1x1_affine.mat",
    )
    assert substitutions[2] == (
        f"sub-ADNI022S0004_ses-M000_{suffix}Warped.nii.gz",
        f"sub-ADNI022S0004_ses-M000_space-MNI152NLin2009cSym_res-1x1x1_{suffix}.nii.gz",
    )
