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
    JHUTract0 = "JHUTracts0"
    JHUTract25 = "JHUTracts25"
    JHUTracts50 = "JHUTracts50"


def _get_resolution_along_axis(label_image_header: Nifti1Header, axis: int) -> str:
    voxels_labels = label_image_header.get_zooms()
    if not 0 <= axis < len(voxels_labels):
        raise ValueError(
            f"The label image has dimension {len(voxels_labels)} and "
            f"axis {axis} is therefore not valid. Please use a value "
            f"between 0 and {len(voxels_labels) - 1}."
        )
    if int(voxels_labels[axis]) == voxels_labels[axis]:
        return str(int(voxels_labels[axis]))
    return str(voxels_labels[axis])


class BaseAtlas:
    """Base class for Atlas handling."""

    __metaclass__ = abc.ABCMeta

    def __init__(self, name: str, roi_filename: str):
        self.name = name
        self.roi_filename = roi_filename
        self.atlas_dir = None
        self.atlas_filename = None

    @property
    @abc.abstractmethod
    def expected_checksum(self) -> str:
        raise NotImplementedError

    @property
    def atlas_folder(self) -> Path:
        return Path(__file__).parent.parent / "resources" / "atlases"

    @property
    def tsv_roi(self) -> Path:
        """Path to the parcellation TSV file.

        The TSV file must contain the `roi_value` and `roi_name` columns:

        roi_value   roi_name
        0   Background
        2001    Precentral_L
        [...]   [...]
        9170    Vermis_10
        """
        return self.atlas_folder / self.roi_filename

    @property
    def spatial_resolution(self) -> str:
        """Return the spatial resolution of the atlas (in format "XxXxX" e.g. 1.5x1.5x1.5)."""
        img_labels = nib.load(self.labels)
        return "x".join(
            _get_resolution_along_axis(img_labels.header, axis=axis)
            for axis in range(3)
        )

    @property
    def labels(self) -> Path:
        """Path to the parcellation in NIfTI format.

        Raises
        ------
        IOError :
            If the checksum of the parcellation found is different from
            the expected checksum.
        """
        from .inputs import compute_sha256_hash

        atlas_labels = self.atlas_dir / self.atlas_filename
        if (checksum := compute_sha256_hash(atlas_labels)) != self.expected_checksum:
            raise IOError(
                f"{atlas_labels} has an SHA256 checksum ({checksum}) "
                f"differing from expected ({self.expected_checksum}), "
                f"file may be corrupted and changed with newer version of FSL."
            )
        return atlas_labels

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
        self.atlas_dir = fsl_dir / "data" / "atlases" / "JHU"


class JHUDTI811mm(FSLAtlas):
    """JHUDTI811mm atlas.

    This atlas contains 48 white matter tract labels that were created by manually
    segmenting a standard-space average of diffusion MRI tensor maps from 81 subjects.

    References
    ----------
    https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/Atlases
    https://www.sciencedirect.com/science/article/abs/pii/S105381190700688X?via%3Dihub
    """

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
    """JHUTracts01mm atlas.

    This atlas contains 20 white matter tract labels that were identified probabilistically
    by averaging the results of deterministic tractography run on 28 subjects.
    Threshold used is 0%.

    References
    ----------
    https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/Atlases
    https://shop.elsevier.com/books/mri-atlas-of-human-white-matter/mori/978-0-444-51741-8
    """

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
    """JHUTracts251mm atlas.

    This atlas contains 20 white matter tract labels that were identified probabilistically
    by averaging the results of deterministic tractography run on 28 subjects.
    Threshold used is 25%.

    References
    ----------
    https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/Atlases
    https://shop.elsevier.com/books/mri-atlas-of-human-white-matter/mori/978-0-444-51741-8
    """

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
    """JHUTracts501mm atlas.

    This atlas contains 20 white matter tract labels that were identified probabilistically
    by averaging the results of deterministic tractography run on 28 subjects.
    Threshold used is 50%.

    References
    ----------
    https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/Atlases
    https://shop.elsevier.com/books/mri-atlas-of-human-white-matter/mori/978-0-444-51741-8
    """

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
    """Local atlases will look for labels in the local 'resources' folder.

    More precisely, the labels and TSV files associated with the atlas
    are located in the folder `<clinica>/resources/atlases/`.
    """

    def __init__(self, name: str, roi_filename: str, atlas_filename: str):
        super().__init__(name, roi_filename)
        self.atlas_filename = atlas_filename
        self.atlas_dir = Path(__file__).parent.parent / "resources" / "atlases"


class AAL2(LocalAtlas):
    """AAL2 atlas.

    Anatomical atlas based on a single subject. It is the updated version of AAL, which is
    probably the most widely used cortical parcellation map in the neuroimaging literature.
    It was built using manual tracing on the spatially normalized single-subject high-resolution
    T1 volume in MNI space. It is composed of 120 regions covering the whole cortex as well as
    the main subcortical structures.

    References
    ----------
    https://www.gin.cnrs.fr/en/tools/aal/
    https://www.sciencedirect.com/science/article/abs/pii/S1053811901909784?via%3Dihub
    """

    def __init__(self):
        super().__init__(
            name="AAL2",
            roi_filename="atlas-AAL2_dseg.tsv",
            atlas_filename="atlas-AAL2_dseg.nii.gz",
        )

    @property
    def expected_checksum(self) -> str:
        return "f6bc698f778a4b383abd3ce355bfd4505c4aa14708e4a7848f8ee928c2b56b37"


class AICHA(LocalAtlas):
    """AICHA atlas.

    Functional atlas based on multiple subjects. It was built using parcellation of group-level
    functional connectivity profiles computed from resting-state fMRI data of 281 healthy subjects.
    It is composed of 384 regions covering the whole cortex as well as the main subcortical structures.

    References
    ----------
    https://www.gin.cnrs.fr/en/tools/aicha/
    https://www.sciencedirect.com/science/article/abs/pii/S0165027015002678?via%3Dihub
    """

    def __init__(self):
        super().__init__(
            name="AICHA",
            roi_filename="atlas-AICHA_dseg.tsv",
            atlas_filename="atlas-AICHA_dseg.nii.gz",
        )

    @property
    def expected_checksum(self) -> str:
        return "cab554d5f546720e60f61f536f82c3d355b31fadb5a4d3ce6a050a606d7ef761"


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
    """Hammers atlas.

    Anatomical atlas based on multiple subjects. It was built using manual tracing on anatomical
    MRI from 30 healthy subjects. The individual subjects parcellations were then registered to MNI
    space to generate a probabilistic atlas as well as a maximum probability map. The latter was
    used in the present work. It is composed of 69 regions covering the whole cortex as well
    as he main subcortical structures.

    References
    ----------
    https://neuro-jena.github.io/cat//index.html#DOWNLOAD
    https://onlinelibrary.wiley.com/doi/epdf/10.1002/hbm.10123
    """

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
    """LPBA40 atlas.

    Anatomical atlas based on multiple subjects. It was built using manual tracing on anatomical
    MRI from 40 healthy subjects. The individual subjects parcellations were then registered to MNI
    space to generate a maximum probability map. It is composed of 56 regions covering the whole
    cortex as well as the main subcortical structures.

    References
    ----------
    https://neuro-jena.github.io/cat//index.html#DOWNLOAD
    https://www.sciencedirect.com/science/article/abs/pii/S1053811907008099?via%3Dihub
    """

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
    """Neuromorphometrics atlas.

    Anatomical atlas based on multiple subjects. It was built using manual tracing on anatomical
    MRI from 30 healthy subjects. The individual subjects parcellations were then registered to
    MNI space to generate a maximum probability map. It is composed of 140 regions covering the
    whole cortex as well as the main subcortical structures. Data were made available for the
    “MICCAI 2012 Grand Challenge and Workshop on Multi-Atlas Labeling”.

    References
    ----------
    https://neuro-jena.github.io/cat//index.html#DOWNLOAD
    http://masiweb.vuse.vanderbilt.edu/workshop2012/index.php/Challenge_Details
    """

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
        return AAL2()
    if atlas_name == AtlasName.AICHA:
        return AICHA()
    if atlas_name == AtlasName.HAMMERS:
        return Hammers()
    if atlas_name == AtlasName.LPBA40:
        return LPBA40()
    if atlas_name == AtlasName.NEUROMORPHOMETRICS:
        return Neuromorphometrics()
    if atlas_name == AtlasName.JHUDTI81:
        return JHUDTI811mm()
    if atlas_name == AtlasName.JHUTract0:
        return JHUTracts01mm()
    if atlas_name == AtlasName.JHUTract25:
        return JHUTracts251mm()
    if atlas_name == AtlasName.JHUTracts50:
        return JHUTracts501mm()
