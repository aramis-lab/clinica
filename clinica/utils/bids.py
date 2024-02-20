from dataclasses import dataclass
from enum import Enum
from os import PathLike
from pathlib import Path
from typing import Dict, Tuple, Union


class Extension(str, Enum):
    """Possible extensions in BIDS file names."""

    NIIGZ = ".nii.gz"
    NII = ".nii"
    JSON = ".json"
    TSV = ".tsv"
    MAT = ".mat"
    BVAL = ".bval"
    BVEC = ".bvec"


class Suffix(str, Enum):
    """Possible suffixes in BIDS file names."""

    DWI = "dwi"
    PET = "pet"
    T1W = "t1w"
    T2W = "t2w"
    FLAIR = "flair"
    AFFINE = "affine"
    PROBABILITY = "probability"
    DEFORMATION = "deformation"
    PHASEDIFF = "phasediff"
    MAGNITUDE1 = "magnitude1"
    BRAINMASK = "brainmask"
    STATISTICS = "statistics"


class BIDSLabel(str):
    """A BIDS label is a short string which does not contain symbols
    used for separating entities in BIDS terminology.
    """

    _min_size = 1
    _max_size = 100
    _forbidden_symbols = ("-", "_", ".")

    def __new__(cls, string):
        if not cls._min_size <= len(string) <= cls._max_size:
            raise ValueError(
                f"A string must be between {cls._min_size} and {cls._max_size} to be a valid BIDS label."
            )
        if any([symbol in string for symbol in cls._forbidden_symbols]):
            raise ValueError(
                f"Provided string '{string}' is not a valid BIDS label because "
                f"it contains at least one of these characters: {cls._forbidden_symbols}."
            )
        instance = super().__new__(cls, string)
        return instance


@dataclass
class BIDSFileName:
    """Class modeling a file name following the BIDS specifications."""

    _subject: BIDSLabel
    _session: BIDSLabel
    _suffix: Suffix
    _extension: Extension
    entities: Dict[BIDSLabel, BIDSLabel]

    @property
    def subject(self) -> str:
        return self._subject

    @subject.setter
    def subject(self, subject: str):
        self._subject = BIDSLabel(subject)

    @property
    def session(self) -> str:
        return self._session

    @session.setter
    def session(self, session: str):
        self._session = BIDSLabel(session)

    @property
    def suffix(self) -> str:
        return self._suffix.value

    @suffix.setter
    def suffix(self, suffix: Union[str, Suffix]):
        self._suffix = Suffix(suffix)

    @property
    def extension(self) -> str:
        return self._extension.value

    @extension.setter
    def extension(self, extension: Union[str, Extension]):
        self._extension = Extension(extension)

    @property
    def sub_ses_id(self) -> str:
        return f"sub-{self.subject}_ses-{self.session}"

    @property
    def name(self) -> str:
        if self.entities:
            txt = "_".join([f"{k}-{v}" for k, v in self.entities.items()])
            return f"{self.sub_ses_id}_{txt}_{self.suffix}{self.extension}"
        return f"{self.sub_ses_id}_{self.suffix}{self.extension}"

    @classmethod
    def from_name(cls, filename: Union[str, PathLike]):
        filename, extension = split_name_from_extension(filename)
        entities, suffix = _tokenize_filename_no_ext(filename)
        subject = entities.pop("sub")
        session = entities.pop("ses")
        return cls(
            BIDSLabel(subject),
            BIDSLabel(session),
            Suffix(suffix),
            Extension(extension),
            {BIDSLabel(k): BIDSLabel(v) for k, v in entities.items()},
        )

    def update_entity(self, entity_name: str, entity_value: str):
        self.entities[BIDSLabel(entity_name)] = BIDSLabel(entity_value)

    def delete_entity(self, entity_name: str):
        entity_name = BIDSLabel(entity_name)
        if entity_name in self.entities:
            self.entities.pop(entity_name)


def _tokenize_filename_no_ext(
    filename_without_extension: str,
) -> Tuple[Dict[str, str], str]:
    if "_" not in filename_without_extension:
        raise ValueError(
            f"BIDS file names have entities separated by '_'. "
            f"You provided {filename_without_extension}."
        )
    tokens = filename_without_extension.split("_")
    if len(tokens) < 3:
        raise ValueError(
            f"A valid BIDS filename should have at least 'sub-XXX_ses-YYY_suffix'. "
            f"You provided {filename_without_extension}."
        )
    suffix = tokens.pop()
    if "-" in suffix:
        raise ValueError(
            f"When tokenizing the filename {filename_without_extension}, the suffix "
            f"found was '{suffix}'. It is invalid because it should not contain a '-' symbol."
        )
    if not all(["-" in token for token in tokens]):
        raise ValueError(
            "The BIDS entities should be key-value pairs separated by a '-' symbol."
            f"The entities found are: {tokens}."
        )
    entities = {k: v for k, v in [s.split("-") for s in tokens]}
    return entities, suffix


def split_name_from_extension(filename: Union[str, PathLike]) -> Tuple[str, str]:
    extension = ""
    filename = Path(filename)
    while "." in filename.name:
        extension = filename.suffix + extension
        filename = Path(filename.stem)
    return filename.name, extension
