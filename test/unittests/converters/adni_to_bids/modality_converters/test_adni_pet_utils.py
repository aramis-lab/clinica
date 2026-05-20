import pytest

from clinica.converters.adni_to_bids._modality import (
    ADNIPETPreprocessingStep,
)
from clinica.converters.adni_to_bids.modality_converters._pet_utils import (
    ADNITracer,
)


@pytest.mark.parametrize(
    "tracer, step, expected",
    [
        (ADNITracer.PIB, ADNIPETPreprocessingStep.STEP0, "ADNI Brain PET: Raw PIB"),
        (ADNITracer.FDG, ADNIPETPreprocessingStep.STEP0, "ADNI Brain PET: Raw FDG"),
        (ADNITracer.FDG, ADNIPETPreprocessingStep.STEP1, "Co-registered Dynamic"),
        (ADNITracer.PIB, ADNIPETPreprocessingStep.STEP2, "PIB Co-registered, Averaged"),
    ],
)
def test_define_pet_processing_step_with_tracer(tracer, step, expected):
    from clinica.converters.adni_to_bids.modality_converters._pet_utils import (
        define_pet_processing_step_with_tracer,
    )

    assert expected == define_pet_processing_step_with_tracer(tracer, step)
