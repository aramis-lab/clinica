import logging

from clinica.converters.aibl_to_bids._converter import _warn_about_missing_files


def test_warn_about_missing_files_reports_only_missing_outputs(tmp_path, caplog):
    existing = tmp_path / "converted.nii.gz"
    existing.touch()
    missing_t1 = tmp_path / "missing-t1.nii.gz"
    missing_pet = tmp_path / "missing-pet.nii.gz"

    with caplog.at_level(logging.WARNING, logger="clinica"):
        _warn_about_missing_files([[existing, missing_t1], [None, missing_pet]])

    assert len(caplog.records) == 1
    warning = caplog.records[0]
    assert warning.name == "clinica"
    assert warning.levelno == logging.WARNING
    assert warning.getMessage().splitlines()[1:] == [str(missing_t1), str(missing_pet)]


def test_warn_about_missing_files_is_quiet_when_outputs_exist(tmp_path, caplog):
    existing = tmp_path / "converted.nii.gz"
    existing.touch()

    with caplog.at_level(logging.WARNING, logger="clinica"):
        _warn_about_missing_files([[existing], [None]])

    assert not caplog.records
