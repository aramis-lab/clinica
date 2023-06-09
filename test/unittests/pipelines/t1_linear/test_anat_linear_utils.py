import pytest


@pytest.mark.parametrize("pipeline_name", ["t1-linear", "flair-linear", "fooo"])
def test_get_substitutions_datasink(pipeline_name):
    from clinica.pipelines.t1_linear.anat_linear_utils import get_substitutions_datasink

    suffix = "T1w" if pipeline_name == "t1-linear" else "FLAIR"
    bids_filename = f"sub-ADNI022S0004_ses-M000_{suffix}"

    a, b = get_substitutions_datasink(bids_filename, pipeline_name)

    assert a == bids_filename
    assert len(b) == 3
    assert b[0] == (
        f"sub-ADNI022S0004_ses-M000_{suffix}Warped_cropped.nii.gz",
        f"sub-ADNI022S0004_ses-M000_{suffix}_space-MNI152NLin2009cSym_desc-Crop_res-1x1x1_{suffix}.nii.gz",
    )
    assert b[1] == (
        f"sub-ADNI022S0004_ses-M000_{suffix}0GenericAffine.mat",
        f"sub-ADNI022S0004_ses-M000_{suffix}_space-MNI152NLin2009cSym_res-1x1x1_affine.mat",
    )
    assert b[2] == (
        f"sub-ADNI022S0004_ses-M000_{suffix}Warped.nii.gz",
        f"sub-ADNI022S0004_ses-M000_{suffix}_space-MNI152NLin2009cSym_res-1x1x1_{suffix}.nii.gz",
    )
