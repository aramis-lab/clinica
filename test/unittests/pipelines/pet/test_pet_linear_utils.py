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


@pytest.mark.parametrize(
    "background_value,object_value,expected_threshold",
    [
        (-0.01, 0.01, 0.0),
        (0.0, 0.01, 0.0),
        (0.01, 0.01, 0.01),
        (0.6, 0.01, 0.6),
        (-1.0, -1.0, 0.0),
    ],
)
def test_compute_clipping_threshold(
    tmp_path,
    background_value: float,
    object_value: float,
    expected_threshold: float,
):
    from clinica.pipelines.pet.linear.utils import _compute_clipping_threshold
    from clinica.utils.testing_utils import build_test_image_cubic_object

    build_test_image_cubic_object(
        shape=(10, 10, 10),
        background_value=background_value,
        object_value=object_value,
        object_size=4,
    ).to_filename(tmp_path / "test.nii.gz")

    assert _compute_clipping_threshold(tmp_path / "test.nii.gz") == pytest.approx(
        expected_threshold
    )
