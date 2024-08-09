from pathlib import Path
from typing import Callable, Union

from clinica.iotools.bids_utils import StudyName

__all__ = [
    "convert",
    "converter_factory",
    "get_converter_name",
]


def convert(
    study: Union[str, StudyName],
    path_to_dataset: Union[str, Path],
    bids_dir: Union[str, Path],
    path_to_clinical: Union[str, Path],
    **kwargs,
):
    converter_factory(study)(path_to_dataset, bids_dir, path_to_clinical, **kwargs)


def get_converter_name(study: Union[str, StudyName]) -> str:
    study = StudyName(study)
    if study == StudyName.ADNI:
        return "Adni2Bids"
    if study == StudyName.AIBL:
        return "Aibl2Bids"
    if study == StudyName.GENFI:
        return "GenfiToBids"
    if study == StudyName.HABS:
        return "HabsToBids"
    if study == StudyName.NIFD:
        return "Nifd2Bids"
    if study == StudyName.OASIS:
        return "Oasis2Bids"
    if study == StudyName.OASIS3:
        return "Oasis3ToBids"
    if study == StudyName.UKB:
        return "UkbToBids"
    if study == StudyName.IXI:
        return "IxiToBids"


def converter_factory(study: Union[str, StudyName]) -> Callable:
    study = StudyName(study)
    if study == StudyName.ADNI:
        from .adni_to_bids import convert
    if study == StudyName.AIBL:
        from .aibl_to_bids import convert
    if study == StudyName.GENFI:
        from .genfi_to_bids import convert
    if study == StudyName.HABS:
        from .habs_to_bids import convert
    if study == StudyName.NIFD:
        from .nifd_to_bids import convert
    if study == StudyName.OASIS:
        from .oasis_to_bids import convert
    if study == StudyName.OASIS3:
        from .oasis3_to_bids import convert
    if study == StudyName.UKB:
        from .ukb_to_bids import convert
    if study == StudyName.IXI:
        from .ixi_to_bids import convert
    return convert
