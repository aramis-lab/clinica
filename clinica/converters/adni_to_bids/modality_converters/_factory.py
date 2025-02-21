from pathlib import Path
from typing import Callable, Iterable, Optional, Union

from clinica.converters.adni_to_bids.adni_utils import (
    ADNIModality,
    ADNIModalityConverter,
)

__all__ = ["converter_factory", "modality_converter_factory"]


ConverterInterface = Callable[
    [Path, Path, Path, Path, Optional[Iterable[str]], bool, int], None
]


def _get_converters_for_modality(
    modality: ADNIModality,
) -> list[ADNIModalityConverter]:
    """Map a user-facing modality from ADNI to an iterable of adni modality converters."""
    if modality == ADNIModality.T1:
        return [ADNIModalityConverter.T1]
    if modality == ADNIModality.PET_FDG:
        return [ADNIModalityConverter.PET_FDG, ADNIModalityConverter.PET_FDG_UNIFORM]
    if modality == ADNIModality.PET_AMYLOID:
        return [ADNIModalityConverter.PET_PIB, ADNIModalityConverter.PET_AV45]
    if modality == ADNIModality.PET_TAU:
        return [ADNIModalityConverter.PET_TAU]
    if modality == ADNIModality.DWI:
        return [ADNIModalityConverter.DWI]
    if modality == ADNIModality.FLAIR:
        return [ADNIModalityConverter.FLAIR]
    if modality == ADNIModality.FMRI:
        return [ADNIModalityConverter.FMRI]
    if modality == ADNIModality.FMAP:
        return [ADNIModalityConverter.FMAP]


def converter_factory(converter: ADNIModalityConverter) -> ConverterInterface:
    """Returns the converter associated with the provided ADNIModalityConverter variant."""
    if converter == ADNIModalityConverter.T1:
        from ._t1 import convert_t1

        return convert_t1
    if converter == ADNIModalityConverter.PET_FDG:
        from ._fdg_pet import convert_fdg_pet

        return convert_fdg_pet
    if converter == ADNIModalityConverter.PET_FDG_UNIFORM:
        from ._fdg_pet import convert_fdg_pet_uniform

        return convert_fdg_pet_uniform
    if converter == ADNIModalityConverter.PET_PIB:
        from ._pib_pet import convert_pib_pet

        return convert_pib_pet
    if converter == ADNIModalityConverter.PET_AV45:
        from ._av45_fbb_pet import convert_av45_fbb_pet

        return convert_av45_fbb_pet
    if converter == ADNIModalityConverter.PET_TAU:
        from ._tau_pet import convert_tau_pet

        return convert_tau_pet
    if converter == ADNIModalityConverter.DWI:
        from ._dwi import convert_dwi

        return convert_dwi
    if converter == ADNIModalityConverter.FLAIR:
        from ._flair import convert_flair

        return convert_flair
    if converter == ADNIModalityConverter.FMRI:
        from ._fmri import convert_fmri

        return convert_fmri
    if converter == ADNIModalityConverter.FMAP:
        from ._fmap import convert_fmap

        return convert_fmap


def modality_converter_factory(
    modality: Union[str, ADNIModality],
) -> Iterable[ConverterInterface]:
    """Returns an iterable of ADNI converters based on the provided modality."""
    return [
        converter_factory(converter)
        for converter in _get_converters_for_modality(ADNIModality(modality))
    ]
