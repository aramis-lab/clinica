import pytest

from clinica.converters.adni_to_bids._modality import (
    ADNIModalityConverter,
    ADNIPETPreprocessingStep,
)


@pytest.mark.parametrize(
    "modality",
    ["PET_FDG", "PET_PIB", "PET_AV45", "PET_TAU", "PET_FBB"],
)
def test_adni_modality_converter_is_pet_true(modality):
    assert ADNIModalityConverter[f"{modality}"].is_pet


@pytest.mark.parametrize("modality", ["T1", "DWI", "FLAIR", "FMRI", "FMAP"])
def test_adni_modality_converter_is_pet_false(modality):
    assert not ADNIModalityConverter[f"{modality}"].is_pet


@pytest.mark.parametrize(
    "modality, expected",
    [
        ("T1", "anat"),
        ("DWI", "dwi"),
        ("FLAIR", "anat"),
        ("FMRI", "func"),
        ("FMAP", "fmap"),
        ("PET_FDG", "pet"),
        ("PET_PIB", "pet"),
        ("PET_AV45", "pet"),
        ("PET_TAU", "pet"),
    ],
)
def test_get_output_path(modality, expected):
    assert ADNIModalityConverter[f"{modality}"].output_folder == expected


@pytest.mark.parametrize(
    "modality",
    ["PET_FDG", "PET_PIB", "PET_AV45", "PET_TAU", "PET_FBB", "T1"],
)
def test_write_json_sidecar_false(modality):
    assert not ADNIModalityConverter[f"{modality}"].json_sidecar


@pytest.mark.parametrize("modality", ["DWI", "FMRI", "FMAP", "FLAIR"])
def test_write_json_sidecar_true(modality):
    assert ADNIModalityConverter[f"{modality}"].json_sidecar


@pytest.mark.parametrize(
    "modality",
    [
        "PET_FDG",
        "PET_PIB",
        "PET_AV45",
        "PET_TAU",
        "PET_FBB",
        "T1",
        "FLAIR",
    ],
)
def test_should_be_centered_true(modality):
    assert ADNIModalityConverter[f"{modality}"].to_center


@pytest.mark.parametrize("modality", ["DWI", "FMRI", "FMAP"])
def test_should_be_centered_false(modality):
    assert not ADNIModalityConverter[f"{modality}"].to_center


@pytest.mark.parametrize(
    "value, expected",
    [
        (0, "ADNI Brain PET: Raw"),
        (1, "Co-registered Dynamic"),
        (2, "Co-registered, Averaged"),
        (3, "Coreg, Avg, Standardized Image and Voxel Size"),
        (4, "Coreg, Avg, Std Img and Vox Siz, Uniform Resolution"),
        (5, "Coreg, Avg, Std Img and Vox Siz, Uniform 6mm Res"),
    ],
)
def test_adni_preprocessing_step_from_value(value, expected):
    assert expected == ADNIPETPreprocessingStep.from_step_value(value).value


@pytest.mark.parametrize("value", [1.2, "truc"])
def test_adni_preprocessing_step_from_value_error(value):
    with pytest.raises(ValueError):
        ADNIPETPreprocessingStep.from_step_value(value)


@pytest.mark.parametrize(
    "modality, expected",
    [
        ("T1", "_T1w"),
        ("DWI", "_dwi"),
        ("FLAIR", "_FLAIR"),
        ("FMRI", "_task-rest_bold"),
        ("FMAP", "_fmap"),
        ("PET_FDG", "_trc-18FFDG_rec-coregavg_pet"),
        ("PET_PIB", "_trc-11CPIB_rec-coregavg_pet"),
        ("PET_TAU", "_trc-18FAV1451_rec-coregavg_pet"),
    ],
)
def test_get_output_filename(modality, expected):
    from clinica.converters.adni_to_bids._modality import _get_output_filename

    assert _get_output_filename(ADNIModalityConverter[f"{modality}"]) == expected


def test_get_output_filename_with_tracer():
    from clinica.converters.adni_to_bids._modality import _get_output_filename
    from clinica.utils.pet import Tracer

    assert (
        _get_output_filename(ADNIModalityConverter.PET_AV45, tracer=Tracer.AV45)
        == "_trc-18FAV45_rec-coregavg_pet"
    )


def test_get_output_filename_with_step():
    from clinica.converters.adni_to_bids._modality import _get_output_filename

    assert (
        _get_output_filename(
            ADNIModalityConverter.PET_AV45,
            pet_preprocessing_step=ADNIPETPreprocessingStep.STEP4_8MM,
        )
        == "_trc-18FAV45_rec-coregiso8_pet"
    )
