from pathlib import Path

import pytest

from clinica.utils.pet import SUVRReferenceRegion


@pytest.mark.parametrize(
    "cropped,suvr_reference_region,expected",
    [
        (
            False,
            SUVRReferenceRegion.PONS,
            "_space-MNI152NLin2009cSym_res-1x1x1_suvr-pons_pet.nii.gz",
        ),
        (
            True,
            SUVRReferenceRegion.CEREBELLUM_PONS,
            "_space-MNI152NLin2009cSym_desc-Crop_res-1x1x1_suvr-cerebellumPons_pet.nii.gz",
        ),
    ],
)
def test_get_pet_bids_components(cropped, suvr_reference_region, expected):
    from clinica.pipelines.pet.linear.utils import _get_pet_bids_components

    assert _get_pet_bids_components(cropped, suvr_reference_region) == expected


def test_rename(tmp_path):
    from clinica.pipelines.pet.linear.utils import _rename

    file_to_rename = tmp_path / "foo.nii.gz"
    file_to_rename.touch()
    source_file = tmp_path / "sub-01_ses-M000_run-01_pet.json"
    source_file.touch()

    output_file = _rename(
        file_to_rename, str(tmp_path / "sub-01_ses-M000_run-01"), "_pet.nii.gz"
    )
    assert file_to_rename.exists()
    assert source_file.exists()
    assert output_file == tmp_path / "sub-01_ses-M000_run-01_pet.nii.gz"
    assert output_file.exists()


def test_rename_into_caps(tmp_path):
    from clinica.pipelines.pet.linear.utils import rename_into_caps

    pet_file_to_rename = tmp_path / "foobarbaz.nii.gz"
    pet_file_to_rename.touch()
    transformation_file_to_rename = tmp_path / "transformation.mat"
    transformation_file_to_rename.touch()
    source_file = tmp_path / "sub-01_ses-M000_run-01_pet.nii.gz"
    source_file.touch()
    a, b, c = rename_into_caps(
        source_file,
        pet_file_to_rename,
        transformation_file_to_rename,
        SUVRReferenceRegion.PONS,
        True,
        output_dir=tmp_path,  # Force the writing to tmp_path instead of current folder...
    )
    assert (
        a
        == tmp_path
        / "sub-01_ses-M000_run-01_space-MNI152NLin2009cSym_res-1x1x1_suvr-pons_pet.nii.gz"
    )
    assert b == tmp_path / "sub-01_ses-M000_run-01_space-T1w_rigid.mat"
    assert c is None
