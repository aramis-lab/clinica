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

T1_VOLUME_ATLASES = [
    "AAL2",
    "AICHA",
    "Hammers",
    "LPBA40",
    "Neuromorphometrics",
]

PET_VOLUME_ATLASES = [
    "AAL2",
    "AICHA",
    "Hammers",
    "LPBA40",
    "Neuromorphometrics",
]

DWI_DTI_ATLASES = [
    "JHUDTI81",
    "JHUTract0",
    "JHUTract25",
]

VOLUME_ATLASES = list(set(T1_VOLUME_ATLASES + PET_VOLUME_ATLASES + DWI_DTI_ATLASES))


class AtlasAbstract:
    """Abstract class for Atlas handling.

    Naming convention for children classes of AtlasAbstract:
    <name_atlas>[<resolution>][<map>]
    """

    __metaclass__ = abc.ABCMeta

    @staticmethod
    @abc.abstractmethod
    def get_name_atlas():
        """Return the name of the atlas (as defined in BIDS/CAPS specifications)."""

    def get_spatial_resolution(self):
        """Return the spatial resolution of the atlas (in format "XxXxX" e.g. 1.5x1.5x1.5)."""
        import nibabel as nib

        img_labels = nib.load(self.get_atlas_labels())
        voxels_labels = img_labels.header.get_zooms()
        # Will display integers without decimals
        if int(voxels_labels[0]) == voxels_labels[0]:
            s_x = str(int(voxels_labels[0]))
        else:
            s_x = str(voxels_labels[0])
        if int(voxels_labels[1]) == voxels_labels[1]:
            s_y = str(int(voxels_labels[1]))
        else:
            s_y = str(voxels_labels[1])
        if int(voxels_labels[2]) == voxels_labels[2]:
            s_z = str(int(voxels_labels[2]))
        else:
            s_z = str(voxels_labels[2])

        return f"{s_x}x{s_y}x{s_z}"

    @staticmethod
    @abc.abstractmethod
    def get_atlas_labels():
        """Return the image with the different labels/ROIs."""

    @staticmethod
    @abc.abstractmethod
    def get_tsv_roi():
        """Return the TSV file containing the ROI (regions of interest) of the atlas."""

    def get_index(self):
        import nibabel as nib
        import numpy as np

        img_labels = nib.load(self.get_atlas_labels())
        img_labels = img_labels.get_fdata(dtype="float32")
        labels = list(set(img_labels.ravel()))
        index_vector = np.zeros(len(labels))
        for index, n in enumerate(labels):
            index_vector[index] = index
        return index_vector


class JHUDTI811mm(AtlasAbstract):
    def __init__(self):
        AtlasAbstract.__init__(self)

    @staticmethod
    def get_name_atlas():
        return "JHUDTI81"

    @staticmethod
    def get_atlas_labels():
        import os

        import nipype.interfaces.fsl as fsl

        from .check_dependency import check_environment_variable
        from .inputs import _sha256

        fsl_dir = check_environment_variable("FSLDIR", "FSL")
        atlas_labels = os.path.join(
            fsl_dir, "data", "atlases", "JHU", "JHU-ICBM-labels-1mm.nii.gz"
        )

        # Adding checksum for updated file with version 6.0.5 of fsl
        fsl_atlas_checksums = {
            "old": "fac584ec75ff2a8631710d3345df96733ed87d9bde3387f5b462f8d22914ed69",
            "new": "3c3f5d2f1250a3df60982acff35a75b99fd549a05d5f8124a63f78221aa0ec16",
        }

        if ["5", "0", "5"] <= fsl.Info.version().split(".") < ["6", "0", "5"]:
            expected_checksum = fsl_atlas_checksums["old"]
        else:
            expected_checksum = fsl_atlas_checksums["new"]

        if _sha256(atlas_labels) != expected_checksum:
            raise IOError(
                f"{atlas_labels} has an SHA256 checksum ({_sha256(atlas_labels)}) "
                f"differing from expected ({expected_checksum}), "
                f"file may be corrupted and changed with newer version of FSL."
            )
        return atlas_labels

    @staticmethod
    def get_tsv_roi():
        from os.path import join, realpath, split

        return join(
            split(realpath(__file__))[0],
            "..",
            "resources",
            "atlases",
            "atlas-JHUDTI81_dseg.tsv",
        )


class JHUTracts01mm(AtlasAbstract):
    def __init__(self):
        AtlasAbstract.__init__(self)

    @staticmethod
    def get_name_atlas():
        return "JHUTracts0"

    @staticmethod
    def get_atlas_labels():
        import os

        from .check_dependency import check_environment_variable
        from .inputs import _sha256

        fsl_dir = check_environment_variable("FSLDIR", "FSL")
        atlas_labels = os.path.join(
            fsl_dir, "data", "atlases", "JHU", "JHU-ICBM-tracts-maxprob-thr0-1mm.nii.gz"
        )
        expected_checksum = (
            "eb1de9413a46b02d2b5c7b77852097c6f42c8a5d55a5dbdef949c2e63b95354e"
        )
        if _sha256(atlas_labels) != expected_checksum:
            raise IOError(
                f"{atlas_labels} has an SHA256 checksum ({_sha256(atlas_labels)}) "
                f"differing from expected ({expected_checksum}), "
                f"file may be corrupted and changed with newer version of FSL."
            )
        return atlas_labels

    @staticmethod
    def get_tsv_roi():
        from os.path import join, realpath, split

        return join(
            split(realpath(__file__))[0],
            "..",
            "resources",
            "atlases",
            "atlas-JHUTract_dseg.tsv",
        )


class JHUTracts251mm(AtlasAbstract):
    def __init__(self):
        AtlasAbstract.__init__(self)

    @staticmethod
    def get_name_atlas():
        return "JHUTracts25"

    @staticmethod
    def get_atlas_labels():
        import os

        from .check_dependency import check_environment_variable
        from .inputs import _sha256

        fsl_dir = check_environment_variable("FSLDIR", "FSL")
        atlas_labels = os.path.join(
            fsl_dir,
            "data",
            "atlases",
            "JHU",
            "JHU-ICBM-tracts-maxprob-thr25-1mm.nii.gz",
        )
        expected_checksum = (
            "7cd85fa2be1918fc83173e9bc0746031fd4c08d70d6c81b7b9224b5d3da6d8a6"
        )
        if _sha256(atlas_labels) != expected_checksum:
            raise IOError(
                f"{atlas_labels} has an SHA256 checksum ({_sha256(atlas_labels)}) "
                f"differing from expected ({expected_checksum}), "
                f"file may be corrupted and changed with newer version of FSL."
            )
        return atlas_labels

    @staticmethod
    def get_tsv_roi():
        from os.path import join, realpath, split

        return join(
            split(realpath(__file__))[0],
            "..",
            "resources",
            "atlases",
            "atlas-JHUTract_dseg.tsv",
        )


class JHUTracts501mm(AtlasAbstract):
    def __init__(self):
        AtlasAbstract.__init__(self)

    @staticmethod
    def get_name_atlas():
        return "JHUTracts50"

    @staticmethod
    def get_atlas_labels():
        import os

        from .check_dependency import check_environment_variable
        from .inputs import _sha256

        fsl_dir = check_environment_variable("FSLDIR", "FSL")
        atlas_labels = os.path.join(
            fsl_dir,
            "data",
            "atlases",
            "JHU",
            "JHU-ICBM-tracts-maxprob-thr50-1mm.nii.gz",
        )
        expected_checksum = (
            "20ff0216d770686838de26393c0bdac38c8104760631a1a2b5f518bc0bbb470a"
        )
        if _sha256(atlas_labels) != expected_checksum:
            raise IOError(
                f"{atlas_labels} has an SHA256 checksum ({_sha256(atlas_labels)}) "
                f"differing from expected ({expected_checksum}), "
                f"file may be corrupted and changed with newer version of FSL."
            )
        return atlas_labels

    @staticmethod
    def get_tsv_roi():
        from os.path import join, realpath, split

        return join(
            split(realpath(__file__))[0],
            "..",
            "resources",
            "atlases",
            "atlas-JHUTract_dseg.tsv",
        )


class AAL2(AtlasAbstract):
    def __init__(self):
        AtlasAbstract.__init__(self)

    @staticmethod
    def get_name_atlas():
        return "AAL2"

    @staticmethod
    def get_atlas_labels():
        from os.path import join, realpath, split

        return join(
            split(realpath(__file__))[0],
            "..",
            "resources",
            "atlases",
            "atlas-AAL2_dseg.nii.gz",
        )

    @staticmethod
    def get_tsv_roi():
        from os.path import join, realpath, split

        return join(
            split(realpath(__file__))[0],
            "..",
            "resources",
            "atlases",
            "atlas-AAL2_dseg.tsv",
        )


class Hammers(AtlasAbstract):
    def __init__(self):
        AtlasAbstract.__init__(self)

    @staticmethod
    def get_name_atlas():
        return "Hammers"

    @staticmethod
    def get_atlas_labels():
        from clinica.utils.inputs import RemoteFileStructure, get_file_from_server

        hammers_parc = RemoteFileStructure(
            filename="atlas-Hammers_dseg.nii.gz",
            url="https://aramislab.paris.inria.fr/files/software/cat12/CAT12-Atlases/",
            checksum="c034a7bce2dcab390a0b72f4e7d04769eb3fe5b990d0e18d89b0ce73339a5376",
        )
        return get_file_from_server(hammers_parc)

    @staticmethod
    def get_tsv_roi():
        from os.path import join, realpath, split

        return join(
            split(realpath(__file__))[0],
            "..",
            "resources",
            "atlases",
            "atlas-Hammers_dseg.tsv",
        )


class LPBA40(AtlasAbstract):
    def __init__(self):
        AtlasAbstract.__init__(self)

    @staticmethod
    def get_name_atlas():
        return "LPBA40"

    @staticmethod
    def get_atlas_labels():
        from clinica.utils.inputs import RemoteFileStructure, get_file_from_server

        lpba40_parc = RemoteFileStructure(
            filename="atlas-LPBA40_dseg.nii.gz",
            url="https://aramislab.paris.inria.fr/files/software/cat12/CAT12-Atlases/",
            checksum="20826b572bbbdbcdbf28bbd3801dc0c2fed28d1e54bc4fd5027e64ccc6d50374",
        )
        return get_file_from_server(lpba40_parc)

    @staticmethod
    def get_tsv_roi():
        from os.path import join, realpath, split

        return join(
            split(realpath(__file__))[0],
            "..",
            "resources",
            "atlases",
            "atlas-LPBA40_dseg.tsv",
        )


class AICHA(AtlasAbstract):
    def __init__(self):
        AtlasAbstract.__init__(self)

    @staticmethod
    def get_name_atlas():
        return "AICHA"

    @staticmethod
    def get_atlas_labels():
        from os.path import join, realpath, split

        return join(
            split(realpath(__file__))[0],
            "..",
            "resources",
            "atlases",
            "atlas-AICHA_dseg.nii.gz",
        )

    @staticmethod
    def get_tsv_roi():
        from os.path import join, realpath, split

        return join(
            split(realpath(__file__))[0],
            "..",
            "resources",
            "atlases",
            "atlas-AICHA_dseg.tsv",
        )


class Neuromorphometrics(AtlasAbstract):
    def __init__(self):
        AtlasAbstract.__init__(self)

    @staticmethod
    def get_name_atlas():
        return "Neuromorphometrics"

    @staticmethod
    def get_atlas_labels():
        from clinica.utils.inputs import RemoteFileStructure, get_file_from_server

        neuromorphometrics_parc = RemoteFileStructure(
            filename="atlas-Neuromorphometrics_dseg.nii.gz",
            url="https://aramislab.paris.inria.fr/files/software/cat12/CAT12-Atlases/",
            checksum="19a50136cd2f8a14357a19ad8a1dc4a2ecb6beb3fc16cb5441f4f2ebaf64a9a5",
        )
        return get_file_from_server(neuromorphometrics_parc)

    @staticmethod
    def get_tsv_roi():
        from os.path import join, realpath, split

        return join(
            split(realpath(__file__))[0],
            "..",
            "resources",
            "atlases",
            "atlas-Neuromorphometrics_dseg.tsv",
        )
