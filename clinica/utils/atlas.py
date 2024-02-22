"""This module contains utilities to handle atlases in Clinica.

An atlas is currently defined by its name and a set of labels in a template space.

This current implementation has some drawbacks:
- Atlas is misleading: it is only a set of labels in a template space
- This implementation can not handle case where parcellation is in several template space
(e.g. `MNI152NLin2009cAsym` and `MNI152NLin6Asym`, see TemplateFlow for example)

Either a refactoring of this module or the use of an external API
(e.g. TemplateFlow - https://www.templateflow.org/) needs to be considered.
"""

import abc
from enum import Enum
from pathlib import Path
from typing import Union

import nibabel as nib
import numpy as np
from nibabel import Nifti1Header


class T1AndPetVolumeAtlasName(str, Enum):
    """Possible names for T1 / PET atlases."""

    AAL2 = "AAL2"
    AICHA = "AICHA"
    HAMMERS = "Hammers"
    LPBA40 = "LPBA40"
    NEUROMORPHOMETRICS = "Neuromorphometrics"


class AtlasName(str, Enum):
    """Possible names for atlases."""

    AAL2 = "AAL2"
    AICHA = "AICHA"
    HAMMERS = "Hammers"
    LPBA40 = "LPBA40"
    NEUROMORPHOMETRICS = "Neuromorphometrics"
    JHUDTI81 = "JHUDTI81"
    JHUTract0 = "JHUTract0"
    JHUTract25 = "JHUTract25"
    JHUTracts50 = "JHUTracts50"


def _get_resolution_along_axis(label_image_header: Nifti1Header, axis: int) -> str:
    voxels_labels = label_image_header.get_zooms()
    if int(voxels_labels[axis]) == voxels_labels[axis]:
        return str(int(voxels_labels[axis]))
    return str(voxels_labels[axis])


class BaseAtlas:
    """Base class for Atlas handling."""

    __metaclass__ = abc.ABCMeta

    def __init__(self, name: str, roi_filename: str):
        self.name = name
        self.roi_filename = roi_filename

    @property
    @abc.abstractmethod
    def expected_checksum(self) -> str:
        raise NotImplementedError

    @property
    def atlas_folder(self) -> Path:
        return Path(__file__).parent.parent / "resources" / "atlases"

    @property
    def tsv_roi(self) -> Path:
        return self.atlas_folder / self.roi_filename

    @property
    def spatial_resolution(self) -> str:
        """Return the spatial resolution of the atlas (in format "XxXxX" e.g. 1.5x1.5x1.5)."""
        img_labels = nib.load(self.get_atlas_labels())
        return "x".join(
            _get_resolution_along_axis(img_labels.header, axis=axis)
            for axis in range(3)
        )

    @property
    @abc.abstractmethod
    def labels(self) -> Path:
        """Return the image with the different labels/ROIs."""
        raise NotImplementedError

    def get_index(self) -> np.ndarray:
        img_labels = nib.load(self.labels)
        img_labels = img_labels.get_fdata(dtype="float32")
        labels = list(set(img_labels.ravel()))
        index_vector = np.zeros(len(labels))
        for index, n in enumerate(labels):
            index_vector[index] = index
        return index_vector


class FSLAtlas(BaseAtlas):
    """FSL atlases look for the labels in the FSL folder (requires FSL)."""

    def __init__(self, name: str, roi_filename: str, atlas_filename: str):
        from .check_dependency import check_environment_variable

        super().__init__(name, roi_filename)
        self.atlas_filename = atlas_filename
        fsl_dir = Path(check_environment_variable("FSLDIR", "FSL"))
        self.fsl_atlas_dir = fsl_dir / "data" / "atlases" / "JHU"

    @property
    def labels(self) -> Path:
        from .inputs import compute_sha256_hash

        atlas_labels = self.fsl_atlas_dir / self.atlas_filename
        if (checksum := compute_sha256_hash(atlas_labels)) != self.expected_checksum:
            raise IOError(
                f"{atlas_labels} has an SHA256 checksum ({checksum}) "
                f"differing from expected ({self.expected_checksum}), "
                f"file may be corrupted and changed with newer version of FSL."
            )
        return atlas_labels


class JHUDTI811mm(FSLAtlas):
    def __init__(self):
        super().__init__(
            name="JHUDTI81",
            roi_filename="atlas-JHUDTI81_dseg.tsv",
            atlas_filename="JHU-ICBM-labels-1mm.nii.gz",
        )

    @property
    def expected_checksum(self) -> str:
        import nipype.interfaces.fsl as fsl

        if ["5", "0", "5"] <= fsl.Info.version().split(".") < ["6", "0", "5"]:
            return "fac584ec75ff2a8631710d3345df96733ed87d9bde3387f5b462f8d22914ed69"
        return "3c3f5d2f1250a3df60982acff35a75b99fd549a05d5f8124a63f78221aa0ec16"


class JHUTracts01mm(FSLAtlas):
    def __init__(self):
        super().__init__(
            name="JHUTracts0",
            roi_filename="atlas-JHUTract_dseg.tsv",
            atlas_filename="JHU-ICBM-tracts-maxprob-thr0-1mm.nii.gz",
        )

    @property
    def expected_checksum(self) -> str:
        return "eb1de9413a46b02d2b5c7b77852097c6f42c8a5d55a5dbdef949c2e63b95354e"


class JHUTracts251mm(FSLAtlas):
    def __init__(self):
        super().__init__(
            name="JHUTracts25",
            roi_filename="atlas-JHUTract_dseg.tsv",
            atlas_filename="JHU-ICBM-tracts-maxprob-thr25-1mm.nii.gz",
        )

    @property
    def expected_checksum(self) -> str:
        return "7cd85fa2be1918fc83173e9bc0746031fd4c08d70d6c81b7b9224b5d3da6d8a6"


class JHUTracts501mm(FSLAtlas):
    def __init__(self):
        super().__init__(
            name="JHUTracts50",
            roi_filename="atlas-JHUTract_dseg.tsv",
            atlas_filename="JHU-ICBM-tracts-maxprob-thr50-1mm.nii.gz",
        )

    @property
    def expected_checksum(self) -> str:
        return "20ff0216d770686838de26393c0bdac38c8104760631a1a2b5f518bc0bbb470a"


class LocalAtlas(BaseAtlas):
    """Local atlases will look for labels in the local 'resources' folder."""

    def __init__(self, name: str, roi_filename: str, atlas_filename: str):
        super().__init__(name, roi_filename)
        self.atlas_filename = atlas_filename

    @property
    def labels(self) -> Path:
        return Path(__file__).parent / "resources" / "atlases" / self.atlas_filename


class AAL2(LocalAtlas):
    def __init__(self):
        super().__init__(
            name="AAL2",
            roi_filename="atlas-AAL2_dseg.tsv",
            atlas_filename="atlas-AAL2_dseg.nii.gz",
        )


class AICHA(LocalAtlas):
    def __init__(self):
        super().__init__(
            name="AICHA",
            roi_filename="atlas-AICHA_dseg.tsv",
            atlas_filename="atlas-AICHA_dseg.nii.gz",
        )


class RemoteAtlas(BaseAtlas):
    """Remote atlases will download the labels from the aramislab server."""

    def __init__(self, name: str, roi_filename: str, atlas_filename: str):
        super().__init__(name, roi_filename)
        self.atlas_filename = atlas_filename

    @property
    def labels(self) -> Path:
        from clinica.utils.inputs import RemoteFileStructure, get_file_from_server

        return get_file_from_server(
            RemoteFileStructure(
                filename=self.atlas_filename,
                url="https://aramislab.paris.inria.fr/files/software/cat12/CAT12-Atlases/",
                checksum=self.expected_checksum,
            )
        )


class Hammers(RemoteAtlas):
    def __init__(self):
        super().__init__(
            name="Hammers",
            roi_filename="atlas-Hammers_dseg.tsv",
            atlas_filename="atlas-Hammers_dseg.nii.gz",
        )

    @property
    def expected_checksum(self) -> str:
        return "c034a7bce2dcab390a0b72f4e7d04769eb3fe5b990d0e18d89b0ce73339a5376"


class LPBA40(RemoteAtlas):
    def __init__(self):
        super().__init__(
            name="LPBA40",
            roi_filename="atlas-LPBA40_dseg.tsv",
            atlas_filename="atlas-LPBA40_dseg.nii.gz",
        )

    @property
    def expected_checksum(self) -> str:
        return "20826b572bbbdbcdbf28bbd3801dc0c2fed28d1e54bc4fd5027e64ccc6d50374"


class Neuromorphometrics(RemoteAtlas):
    def __init__(self):
        super().__init__(
            name="Neuromorphometrics",
            roi_filename="atlas-Neuromorphometrics_dseg.tsv",
            atlas_filename="atlas-Neuromorphometrics_dseg.nii.gz",
        )

    @property
    def expected_checksum(self) -> str:
        return "19a50136cd2f8a14357a19ad8a1dc4a2ecb6beb3fc16cb5441f4f2ebaf64a9a5"


def atlas_factory(atlas_name: Union[str, AtlasName, BaseAtlas]) -> BaseAtlas:
    """Factory method for atlases.

    Parameters
    ----------
    atlas_name : str or AtlasName or atlas instance
        If an atlas instance, the instance is returned.
        If a string, the corresponding atlas will be returned.

    Returns
    -------
    BaseAtlas :
        The atlas class corresponding to the provided name.
    """
    if isinstance(atlas_name, BaseAtlas):
        return atlas_name
    if isinstance(atlas_name, str):
        atlas_name = AtlasName(atlas_name)
    if atlas_name == AtlasName.AAL2:
        return AAL2
    if atlas_name == AtlasName.AICHA:
        return AICHA
    if atlas_name == AtlasName.HAMMERS:
        return Hammers
    if atlas_name == AtlasName.LPBA40:
        return LPBA40
    if atlas_name == AtlasName.NEUROMORPHOMETRICS:
        return Neuromorphometrics
    if atlas_name == AtlasName.JHUDTI81:
        return JHUDTI811mm
    if atlas_name == AtlasName.JHUTract0:
        return JHUTracts01mm
    if atlas_name == AtlasName.JHUTract25:
        return JHUTracts251mm
    if atlas_name == AtlasName.JHUTracts50:
        return JHUTracts501mm
