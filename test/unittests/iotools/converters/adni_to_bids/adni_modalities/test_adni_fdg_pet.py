import pytest


@pytest.mark.parametrize("step_value,expected", [(2, "fdg"), (4, "fdg_uniform")])
def test_get_modality_from_adni_preprocessing_step(step_value, expected):
    from clinica.iotools.converters.adni_to_bids.adni_modalities.adni_fdg_pet import (
        ADNIPreprocessingStep,
        _get_modality_from_adni_preprocessing_step,
    )

    assert (
        _get_modality_from_adni_preprocessing_step(
            ADNIPreprocessingStep.from_step_value(step_value)
        )
        == expected
    )


@pytest.mark.parametrize("step_value", [0, 1, 3, 5])
def test_get_modality_from_adni_preprocessing_step_error(step_value):
    from clinica.iotools.converters.adni_to_bids.adni_modalities.adni_fdg_pet import (
        ADNIPreprocessingStep,
        _get_modality_from_adni_preprocessing_step,
    )

    with pytest.raises(
        ValueError,
        match="The ADNI preprocessing step",
    ):
        _get_modality_from_adni_preprocessing_step(
            ADNIPreprocessingStep.from_step_value(step_value)
        )
