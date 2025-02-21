import re
from abc import ABC, abstractmethod
from collections import UserString
from enum import Enum
from typing import Type

__all__ = [
    "StudyName",
    "BIDSSubjectID",
    "bids_id_factory",
    "ADNIBIDSSubjectID",
    "NIFDBIDSSubjectID",
    "AIBLBIDSSubjectID",
    "UKBBIDSSubjectID",
    "GENFIBIDSSubjectID",
    "OASISBIDSSubjectID",
    "OASIS3BIDSSubjectID",
    "HABSBIDSSubjectID",
    "IXIBIDSSubjectID",
]


class StudyName(str, Enum):
    """Studies supported by the converters of Clinica."""

    ADNI = "ADNI"
    AIBL = "AIBL"
    GENFI = "GENFI"
    HABS = "HABS"
    NIFD = "NIFD"
    OASIS = "OASIS"
    OASIS3 = "OASIS3"
    UKB = "UKB"
    IXI = "IXI"


class BIDSSubjectID(ABC, UserString):
    """This is the interface that BIDS subject IDs have to implement."""

    def __init__(self, value: str):
        instance = super().__init__(self.validate(value))
        return instance

    @abstractmethod
    def validate(self, value: str) -> str:
        raise NotImplementedError

    @classmethod
    @abstractmethod
    def from_original_study_id(cls, study_id: str) -> str:
        raise NotImplementedError

    @abstractmethod
    def to_original_study_id(self) -> str:
        raise NotImplementedError


def bids_id_factory(study: StudyName) -> Type[BIDSSubjectID]:
    if study == StudyName.ADNI:
        return ADNIBIDSSubjectID
    if study == StudyName.NIFD:
        return NIFDBIDSSubjectID
    if study == StudyName.AIBL:
        return AIBLBIDSSubjectID
    if study == StudyName.UKB:
        return UKBBIDSSubjectID
    if study == StudyName.GENFI:
        return GENFIBIDSSubjectID
    if study == StudyName.OASIS:
        return OASISBIDSSubjectID
    if study == StudyName.OASIS3:
        return OASIS3BIDSSubjectID
    if study == StudyName.HABS:
        return HABSBIDSSubjectID
    if study == StudyName.IXI:
        return IXIBIDSSubjectID


class ADNIBIDSSubjectID(BIDSSubjectID):
    """Implementation for ADNI of the BIDSSubjectIDClass, allowing to go from the source id XXX_S_XXXX
    to a bids id sub-ADNIXXXSXXX and reciprocally."""

    def validate(self, value: str) -> str:
        if re.fullmatch(r"sub-ADNI\d{3}S\d{4}", value):
            return value
        raise ValueError(
            f"BIDS ADNI subject ID {value} is not properly formatted. "
            "Expecting a 'sub-ADNIXXXSXXXX' format."
        )

    @classmethod
    def from_original_study_id(cls, study_id: str) -> str:
        if re.fullmatch(r"\d{3}_S_\d{4}", study_id):
            return "sub-ADNI" + study_id.replace("_", "")
        raise ValueError(
            f"Raw ADNI subject ID {study_id} is not properly formatted. "
            "Expecting a 'XXX_S_XXXX' format."
        )

    def to_original_study_id(self) -> str:
        return "_S_".join(self.split("ADNI")[1].split("S"))


class NIFDBIDSSubjectID(BIDSSubjectID):
    """Implementation for NIFD of the BIDSSubjectIDClass, allowing to go from the source id X_S_XXXX
    to a bids id sub-NIFDXSXXX and reciprocally."""

    def validate(self, value: str) -> str:
        if re.fullmatch(r"sub-NIFD\dS\d{4}", value):
            return value
        raise ValueError(
            f"BIDS NIFD subject ID {value} is not properly formatted. "
            "Expecting a 'sub-NIFDXSXXXX' format."
        )

    @classmethod
    def from_original_study_id(cls, study_id: str) -> str:
        if re.fullmatch(r"\d_S_\d{4}", study_id):
            return "sub-NIFD" + study_id.replace("_", "")
        raise ValueError(
            f"Raw NIFD subject ID {study_id} is not properly formatted. "
            "Expecting a 'X_S_XXXX' format."
        )

    def to_original_study_id(self) -> str:
        return "_S_".join(self.split("NIFD")[1].split("S"))


class AIBLBIDSSubjectID(BIDSSubjectID):
    """Implementation for AIBL of the BIDSSubjectIDClass, allowing to go from the source id Y
    to a bids id sub-ADNIY and reciprocally."""

    def validate(self, value: str) -> str:
        if re.fullmatch(r"sub-AIBL\d*", value):
            return value
        raise ValueError(
            f"BIDS AIBL subject ID {value} is not properly formatted. "
            "Expecting a 'sub-AIBLY' format."
        )

    @classmethod
    def from_original_study_id(cls, study_id: str) -> str:
        if re.fullmatch(r"\d*", study_id):
            return "sub-AIBL" + study_id
        raise ValueError(
            f"Raw AIBL subject ID {study_id} is not properly formatted. "
            "Expecting a 'Y' format where Y is a combination of digits."
        )

    def to_original_study_id(self) -> str:
        return self.split("AIBL")[1]


class UKBBIDSSubjectID(BIDSSubjectID):
    """Implementation for UKB of the BIDSSubjectIDClass, allowing to go from the source id Y
    to a bids id sub-ADNIY and reciprocally."""

    def validate(self, value: str) -> str:
        if re.fullmatch(r"sub-UKB\d*", value):
            return value
        raise ValueError(
            f"BIDS UKB subject ID {value} is not properly formatted. "
            "Expecting a 'sub-UKBY' format."
        )

    @classmethod
    def from_original_study_id(cls, study_id: str) -> str:
        if re.fullmatch(r"\d*", study_id):
            return "sub-UKB" + study_id
        raise ValueError(
            f"Raw UKB subject ID {study_id} is not properly formatted. "
            "Expecting a 'Y' format where Y is a combination of digits."
        )

    def to_original_study_id(self) -> str:
        return self.split("UKB")[1]


class GENFIBIDSSubjectID(BIDSSubjectID):
    """Implementation for GENFI of the BIDSSubjectIDClass, allowing to go from the source id Y
    to a bids id sub-Y and reciprocally."""

    def validate(self, value: str) -> str:
        if re.fullmatch(r"sub-\w*", value):
            return value
        raise ValueError(
            f"BIDS GENFI subject ID {value} is not properly formatted. "
            "Expecting a 'sub-Y' format."
        )

    @classmethod
    def from_original_study_id(cls, study_id: str) -> str:
        if re.fullmatch(r"\w*", study_id):
            return "sub-" + study_id
        raise ValueError(
            f"Raw GENFI subject ID {study_id} is not properly formatted. "
            "Expecting a 'Y' format where Y is a combination of letters and digits."
        )

    def to_original_study_id(self) -> str:
        return self.split("-")[1]


class OASISBIDSSubjectID(BIDSSubjectID):
    """Implementation for OASIS1 of the BIDSSubjectIDClass, allowing to go from the source id OAS1_XXXX_MR1/2
    to a bids id sub-OASIS1XXXX and reciprocally."""

    def validate(self, value: str) -> str:
        if re.fullmatch(r"sub-OASIS1\d{4}", value):
            return value
        raise ValueError(
            f"BIDS OASIS1 subject ID {value} is not properly formatted. "
            "Expecting a 'sub-OASIS1XXXX' format."
        )

    @classmethod
    def from_original_study_id(cls, study_id: str) -> str:
        if re.fullmatch(r"OAS1_\d{4}_MR\d", study_id):
            return "sub-OASIS1" + study_id.split("_")[1]
        raise ValueError(
            f"Raw OASIS1 subject ID {study_id} is not properly formatted. "
            "Expecting a 'OAS1_XXXX_MR1/2' format."
        )

    def to_original_study_id(self) -> str:
        return f"OAS1_{self.split('OASIS1')[1]}_MR1"


class OASIS3BIDSSubjectID(BIDSSubjectID):
    """Implementation for OASIS3 of the BIDSSubjectIDClass, allowing to go from the source id XXXX
    to a bids id sub-OAS3XXXX and reciprocally."""

    def validate(self, value: str) -> str:
        if re.fullmatch(r"sub-OAS3\d{4}", value):
            return value
        raise ValueError(
            f"BIDS OASIS3 subject ID {value} is not properly formatted. "
            "Expecting a 'sub-OAS3XXXX' format."
        )

    @classmethod
    def from_original_study_id(cls, study_id: str) -> str:
        if re.fullmatch(r"OAS3\d{4}", study_id):
            return "sub-" + study_id
        raise ValueError(
            f"Raw OASIS3 subject ID {study_id} is not properly formatted. "
            "Expecting a 'OAS3XXXX' format."
        )

    def to_original_study_id(self) -> str:
        return self.split("-")[1]


class HABSBIDSSubjectID(BIDSSubjectID):
    """Implementation for HABS of the BIDSSubjectIDClass, allowing to go from the source id P_Y
    to a bids id sub-HABSY and reciprocally."""

    def validate(self, value: str) -> str:
        if re.fullmatch(r"sub-HABS\w*", value):
            return value
        raise ValueError(
            f"BIDS HABS subject ID {value} is not properly formatted. "
            "Expecting a 'sub-HABSY' format."
        )

    @classmethod
    def from_original_study_id(cls, study_id: str) -> str:
        if re.fullmatch(r"P_\w*", study_id):
            return study_id.replace("P_", "sub-HABS")
        raise ValueError(
            f"Raw HABS subject ID {study_id} is not properly formatted. "
            "Expecting a 'P_Y' format."
        )

    def to_original_study_id(self) -> str:
        return str(self).replace("sub-HABS", "P_")


class IXIBIDSSubjectID(BIDSSubjectID):
    """Implementation for IXI of the BIDSSubjectIDClass, allowing to go from the source id IXI###
    to a bids id sub-IXI### and reciprocally."""

    def validate(self, value: str) -> str:
        if re.fullmatch(r"sub-IXI\d{3}", value):
            return value
        raise ValueError(
            f"BIDS IXI subject ID {value} is not properly formatted. "
            "Expecting a 'sub-IXIXXX' format."
        )

    @classmethod
    def from_original_study_id(cls, study_id: str) -> str:
        if re.fullmatch(r"IXI\d{3}", study_id):
            return f"sub-{study_id}"
        raise ValueError(
            f"Raw IXI subject ID {study_id} is not properly formatted. "
            "Expecting a 'Y' format."
        )

    def to_original_study_id(self) -> str:
        return str(self).replace("sub-", "")
