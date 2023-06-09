import pytest


@pytest.mark.parametrize("pipeline_name", ["t1-linear", "flair-linear", "fooo"])
def test_get_substitutions_datasink(pipeline_name):
    from clinica.pipelines.t1_linear.anat_linear_utils import get_substitutions_datasink

    suffix = "T1w" if pipeline_name == "t1-linear" else "FLAIR"
    bids_image_id = f"sub-ADNI022S0004_ses-M000_{suffix}"

    substitutions = get_substitutions_datasink(bids_image_id, pipeline_name)

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
