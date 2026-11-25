import warnings

import pytest

from clinica.converters.aibl_to_bids._converter import _warn_about_missing_files


def test_warn_about_missing_files_reports_only_missing_outputs(tmp_path):
    existing = tmp_path / "converted.nii.gz"
    existing.touch()
    missing_t1 = tmp_path / "missing-t1.nii.gz"
    missing_pet = tmp_path / "missing-pet.nii.gz"

    with pytest.warns(UserWarning) as record:
        _warn_about_missing_files([[existing, missing_t1], [None, missing_pet]])

    assert len(record) == 1
    assert str(record[0].message).splitlines()[1:] == [
        str(missing_t1),
        str(missing_pet),
    ]


def test_warn_about_missing_files_is_quiet_when_outputs_exist(tmp_path):
    existing = tmp_path / "converted.nii.gz"
    existing.touch()

    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter("always")
        _warn_about_missing_files([[existing], [None]])

    assert not record
