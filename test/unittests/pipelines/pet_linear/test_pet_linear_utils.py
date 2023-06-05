from pathlib import Path

import pytest


@pytest.mark.parametrize(
    "cropped,suvr_reference_region,expected_suffix",
    [
        (False, "foo", "_space-MNI152NLin2009cSym_res-1x1x1_suvr-foo_pet.nii.gz"),
        (
            True,
            "bar",
            "_space-MNI152NLin2009cSym_desc-Crop_res-1x1x1_suvr-bar_pet.nii.gz",
        ),
    ],
)
def test_get_pet_suffix(cropped, suvr_reference_region, expected_suffix):
    from clinica.pipelines.pet_linear.pet_linear_utils import _get_pet_suffix

    assert _get_pet_suffix(cropped, suvr_reference_region) == expected_suffix


def test_rename(tmp_path):
    from clinica.pipelines.pet_linear.pet_linear_utils import _rename

    file_to_rename = tmp_path / "foo.nii.gz"
    file_to_rename.touch()
    source_file = tmp_path / "sub-01_ses-M000_run-01_pet.json"
    source_file.touch()

    output_file = _rename(
        str(file_to_rename), "sub-01_ses-M000_run-01_pet", "_pet.nii.gz"
    )
    assert file_to_rename.exists()
    assert source_file.exists()
    assert Path(output_file).name == "sub-01_ses-M000_run-01_pet_pet.nii.gz"
    assert Path(output_file).exists()


def test_rename_into_caps(tmp_path):
    from clinica.pipelines.pet_linear.pet_linear_utils import rename_into_caps

    pet_file_to_rename = tmp_path / "foobarbaz.nii.gz"
    pet_file_to_rename.touch()
    transformation_file_to_rename = tmp_path / "transformation.mat"
    transformation_file_to_rename.touch()
    source_file = tmp_path / "sub-01_ses-M000_run-01_pet.nii.gz"
    source_file.touch()
    a, b, c = rename_into_caps(
        str(source_file),
        str(pet_file_to_rename),
        str(transformation_file_to_rename),
        "suvrfoo",
        True,
    )
    assert (
        Path(a).name
        == "sub-01_ses-M000_run-01_pet_space-MNI152NLin2009cSym_res-1x1x1_suvr-suvrfoo_pet.nii.gz"
    )
    assert Path(b).name == "sub-01_ses-M000_run-01_pet_space-T1w_rigid.mat"
    assert c is None
