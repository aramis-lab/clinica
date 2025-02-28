import pytest

from clinica.utils.pet import Tracer


@pytest.mark.parametrize(
    "filename, expected",
    [
        ("sub-1_ses-1_trc-18FFDG_pet.nii.gz", Tracer.FDG),
        ("sub-1_ses-1_acq-18FFDG_pet.nii.gz", Tracer.FDG),
        ("sub-1_ses-1_trc-18FFDG_acq-18FFDG_pet.nii.gz", Tracer.FDG),
        ("sub-1_ses-1_trc-11CPIB_pet.nii.gz", Tracer.PIB),
        ("sub-1_ses-1_trc-18FFBB_pet.nii.gz", Tracer.FBB),
        ("sub-1_ses-1_trc-18FFMM_pet.nii.gz", Tracer.FMM),
        ("sub-1_ses-1_trc-18FAV1451_pet.nii.gz", Tracer.AV1451),
        ("sub-1_ses-1_trc-18FAV45_pet.nii.gz", Tracer.AV45),
    ],
)
def test_get_pet_tracer_from_filename_success(filename, expected):
    from clinica.utils.pet import get_pet_tracer_from_filename

    assert expected == get_pet_tracer_from_filename(filename)


@pytest.mark.parametrize(
    "filename",
    [
        "sub-1_ses-1_trc18FFDG_pet.nii.gz",
        "sub-1_ses-1_trc-///_pet.nii.gz",
        "sub-1_ses-1_trc-foo_pet.nii.gz",
        "sub-1_ses-1_pet.nii.gz",
    ],
)
def test_get_pet_tracer_from_filename_error(filename):
    from clinica.utils.pet import get_pet_tracer_from_filename

    with pytest.raises(ValueError):
        get_pet_tracer_from_filename(filename)
