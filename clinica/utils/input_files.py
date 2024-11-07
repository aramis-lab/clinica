"""This module contains dictionaries used in inputs.py::clinica_{file|group}_reader().

These dictionaries describe files to grab.
"""

import functools
from collections.abc import Iterable
from dataclasses import dataclass
from enum import Enum, auto
from pathlib import Path
from typing import Callable, Optional, Union

from clinica.utils.dwi import DTIBasedMeasure
from clinica.utils.image import HemiSphere
from clinica.utils.pet import ReconstructionMethod, SUVRReferenceRegion, Tracer

from .spm import get_spm_tissue_from_index

__all__ = [
    "DWIFileType",
    "Parcellation",
    "QueryPattern",
    "QueryPatternName",
    "query_pattern_factory",
    "get_dwi_file",
    "get_dwi_preprocessed_file",
    "get_dwi_fmap_phasediff_file",
    "get_dwi_fmap_magnitude1_file",
    "get_t1w_linear",
    "get_t1_freesurfer_white_matter_surface",
    "get_t1_freesurfer_longitudinal_white_matter_surface",
    "get_t1_freesurfer_segmentation",
    "get_t1_freesurfer_parcellation",
    "get_t1_freesurfer_template",
    "get_t1_freesurfer_longitudinal_parcellation",
    "get_t1_volume_native_tpm",
    "get_t1_volume_dartel_input_tissue",
]


@dataclass
class QueryPattern:
    """Represents a pattern to be used by the clinica_file_reader to query some specific files.

    Attributes
    ----------
    pattern : str
        The actual pattern string to be used to match file names.

    description : str
        A plain text description of the files the pattern matches.

    needed_pipeline : str
        The pipelines that should have been run in order to have the requested files.
        TODO: Improve this to be an iterable of PipelineName objects.
              The difficult part is that some pattern have combinations with AND and OR.
    """

    pattern: str
    description: str
    needed_pipeline: str

    def __post_init__(self):
        if len(self.pattern) == 0:
            raise ValueError("Pattern cannot be empty.")
        if self.pattern[0] == "/":
            raise ValueError(
                "pattern argument cannot start with char: / (does not work in os.path.join function). "
                "If you want to indicate the exact name of the file, use the format "
                "directory_name/filename.extension or filename.extension in the pattern argument."
            )

    def to_dict(self) -> dict:
        return {
            "pattern": self.pattern,
            "description": self.description,
            "needed_pipeline": self.needed_pipeline,
        }


class QueryPatternName(Enum):
    """The different names for usual pattern in Clinica.

    T1W : Get T1W MRI in BIDS
    T2W : Get T2W FLAIR MRI in BIDS
    T1_FS_WM : GET Freesurfer segmentation of white matter
    T1_FS_BRAIN :  Get Freesurfer extracted brain from T1w MRI
    T1_FS_ORIG_NU : Get Freesurfer intensity normalized volume after correction for non-uniformity
    T1_FS_WM_SURF : Get white matter border surface files from the Freesurfer output
    T1_FS_LONG_SURF : Get white matter border surface files from the Freesurfer longitudinal output
    """

    T1W = auto()
    T1W_LINEAR = auto()
    T1W_TO_MNI_TRANSFORM = auto()
    T2W = auto()
    T1_FREESURFER_WHITE_MATTER = auto()
    T1_FREESURFER_BRAIN = auto()
    T1_FREESURFER_ORIG_NU = auto()
    T1_FREESURFER_LONG_ORIG_NU = auto()
    T1_FREESURFER_WHITE_MATTER_SURFACE = auto()
    T1_FREESURFER_LONG_SURFACE = auto()
    T1_FREESURFER_PARCELLATION = auto()
    T1_FREESURFER_LONG_PARCELLATION = auto()
    T1_FREESURFER_SEGMENTATION = auto()
    T1_FREESURFER_TEMPLATE = auto()
    T1_VOLUME_NATIVE_TPM = auto()
    T1_VOLUME_DARTEL_INPUT_TISSUE = auto()
    T1_VOLUME_DEFORMATION_TO_TEMPLATE = auto()
    T1_VOLUME_GROUP_TEMPLATE = auto()
    T1_VOLUME_ITERATION_GROUP_TEMPLATE = auto()
    T1_VOLUME_TEMPLATE_TPM_IN_MNI = auto()
    DWI = auto()
    DWI_PREPROC = auto()
    DWI_PREPROC_BRAINMASK = auto()
    DWI_FMAP_PHASEDIFF = auto()
    DWI_FMAP_MAGNITUDE1 = auto()
    DWI_DTI = auto()
    PET_NII = auto()
    PET_LINEAR_NII = auto()
    PET_VOLUME_NORMALIZED_SUVR = auto()


class Parcellation(str, Enum):
    """The possible atlas names used for deriving parcellations and segmentations."""

    DESIKAN = "Desikan"
    DESTRIEUX = "Destrieux"


class DWIFileType(str, Enum):
    NII = "nii"
    JSON = "json"
    BVEC = "bvec"
    BVAL = "bval"


QueryPatternBuilderInterface = Callable[..., QueryPattern]


def query_pattern_factory(
    name: Union[str, QueryPatternName],
) -> QueryPatternBuilderInterface:
    """Return the query pattern builder corresponding to the provided name.

    Parameters
    ----------
    name : str or QueryPatternName
        The name of the desired pattern.

    Returns
    -------
    QueryPatternBuilderInterface :
        The desired query pattern builder.
    """
    name = QueryPatternName(name)
    if name == QueryPatternName.T1W:
        return get_t1w_mri
    if name == QueryPatternName.T2W:
        return get_t2w_mri
    if name == QueryPatternName.T1_FREESURFER_WHITE_MATTER:
        return get_t1_freesurfer_segmentation_white_matter
    if name == QueryPatternName.T1_FREESURFER_BRAIN:
        return get_t1_freesurfer_extracted_brain
    if name == QueryPatternName.T1_FREESURFER_ORIG_NU:
        return get_t1_freesurfer_intensity_normalized_volume_after_nu
    if name == QueryPatternName.T1_FREESURFER_LONG_ORIG_NU:
        return get_t1_freesurfer_longitudinal_intensity_normalized_volume_after_nu
    if name == QueryPatternName.T1_FREESURFER_WHITE_MATTER_SURFACE:
        return get_t1_freesurfer_white_matter_surface
    if name == QueryPatternName.T1_FREESURFER_LONG_SURFACE:
        return get_t1_freesurfer_longitudinal_white_matter_surface
    if name == QueryPatternName.T1_VOLUME_NATIVE_TPM:
        return get_t1_volume_native_tpm
    if name == QueryPatternName.T1_VOLUME_DARTEL_INPUT_TISSUE:
        return get_t1_volume_dartel_input_tissue
    if name == QueryPatternName.T1_VOLUME_DEFORMATION_TO_TEMPLATE:
        return get_t1_volume_deformation_to_template
    if name == QueryPatternName.T1_VOLUME_GROUP_TEMPLATE:
        return get_t1_volume_group_template
    if name == QueryPatternName.T1_VOLUME_ITERATION_GROUP_TEMPLATE:
        return get_t1_volume_i_th_iteration_group_template
    if name == QueryPatternName.T1_VOLUME_TEMPLATE_TPM_IN_MNI:
        return get_t1_volume_template_tpm_in_mni
    if name == QueryPatternName.T1W_LINEAR:
        return get_t1w_linear
    if name == QueryPatternName.T1W_TO_MNI_TRANSFORM:
        return get_t1w_to_mni_transform
    if name == QueryPatternName.T1_FREESURFER_PARCELLATION:
        return get_t1_freesurfer_parcellation
    if name == QueryPatternName.T1_FREESURFER_LONG_PARCELLATION:
        return get_t1_freesurfer_longitudinal_parcellation
    if name == QueryPatternName.T1_FREESURFER_SEGMENTATION:
        return get_t1_freesurfer_segmentation
    if name == QueryPatternName.T1_FREESURFER_TEMPLATE:
        return get_t1_freesurfer_template
    if name == QueryPatternName.DWI:
        return get_dwi_file
    if name == QueryPatternName.DWI_PREPROC:
        return get_dwi_preprocessed_file
    if name == QueryPatternName.DWI_PREPROC_BRAINMASK:
        return get_dwi_preprocessed_brainmask
    if name == QueryPatternName.DWI_FMAP_PHASEDIFF:
        return get_dwi_fmap_phasediff_file
    if name == QueryPatternName.DWI_FMAP_MAGNITUDE1:
        return get_dwi_fmap_magnitude1_file
    if name == QueryPatternName.DWI_DTI:
        return get_dwi_dti
    if name == QueryPatternName.PET_NII:
        return get_pet_nifti
    if name == QueryPatternName.PET_LINEAR_NII:
        return get_pet_linear_nifti
    if name == QueryPatternName.PET_VOLUME_NORMALIZED_SUVR:
        return get_pet_volume_normalized_suvr


def get_t1w_mri(*args, **kwargs) -> QueryPattern:
    """Get T1W MRI in BIDS."""
    return QueryPattern("sub-*_ses-*_t1w.nii*", "T1w MRI", "")


def get_t2w_mri(*args, **kwargs) -> QueryPattern:
    """Get T2W FLAIR MRI in BIDS."""
    return QueryPattern("sub-*_ses-*_flair.nii*", "FLAIR T2w MRI", "")


def get_t1_freesurfer_segmentation_white_matter(*args, **kwargs) -> QueryPattern:
    """GET Freesurfer segmentation of white matter."""
    return QueryPattern(
        "t1/freesurfer_cross_sectional/sub-*_ses-*/mri/wm.seg.mgz",
        "segmentation of white matter (mri/wm.seg.mgz).",
        "t1-freesurfer",
    )


def get_t1_freesurfer_extracted_brain(*args, **kwargs) -> QueryPattern:
    return QueryPattern(
        "t1/freesurfer_cross_sectional/sub-*_ses-*/mri/brain.mgz",
        "extracted brain from T1w MRI (mri/brain.mgz).",
        "t1-freesurfer",
    )


def get_t1_freesurfer_intensity_normalized_volume_after_nu(
    *args, **kwargs
) -> QueryPattern:
    return QueryPattern(
        "t1/freesurfer_cross_sectional/sub-*_ses-*/mri/orig_nu.mgz",
        (
            "intensity normalized volume generated after correction for "
            "non-uniformity in FreeSurfer (mri/orig_nu.mgz)."
        ),
        "t1-freesurfer",
    )


def get_t1_freesurfer_longitudinal_intensity_normalized_volume_after_nu(
    *args, **kwargs
) -> QueryPattern:
    return QueryPattern(
        "t1/long-*/freesurfer_longitudinal/sub-*_ses-*.long.sub-*_*/mri/orig_nu.mgz",
        (
            "intensity normalized volume generated after correction for "
            "non-uniformity in FreeSurfer (orig_nu.mgz) in longitudinal"
        ),
        "t1-freesurfer and t1-freesurfer longitudinal",
    )


def get_t1w_to_mni_transform(*args, **kwargs) -> QueryPattern:
    return QueryPattern(
        "*space-MNI152NLin2009cSym_res-1x1x1_affine.mat",
        "Transformation matrix from T1W image to MNI space using t1-linear pipeline",
        "t1-linear",
    )


def get_dwi_file(filetype: Union[str, DWIFileType]) -> QueryPattern:
    """Return the query to get DWI files (nii, json, bvec, bval)."""
    filetype = DWIFileType(filetype)
    return QueryPattern(
        f"dwi/sub-*_ses-*_dwi.{filetype.value}*", f"DWI {filetype.value} files.", ""
    )


def get_dwi_preprocessed_file(filetype: Union[str, DWIFileType]) -> QueryPattern:
    filetype = DWIFileType(filetype)
    return QueryPattern(
        f"dwi/preprocessing/sub-*_ses-*_space-*_desc-preproc_dwi.{filetype.value}*",
        f"preprocessed {filetype.value} files",
        "dwi-preprocessing-using-t1 or dwi-preprocessing-using-fieldmap",
    )


def get_dwi_preprocessed_brainmask(*args, **kwargs) -> QueryPattern:
    return QueryPattern(
        "dwi/preprocessing/sub-*_ses-*_space-*_brainmask.nii*",
        "b0 brainmask",
        "dwi-preprocessing-using-t1 or dwi-preprocessing-using-fieldmap",
    )


def get_dwi_fmap_phasediff_file(filetype: Union[str, DWIFileType]) -> QueryPattern:
    filetype = DWIFileType(filetype)
    return QueryPattern(
        f"fmap/sub-*_ses-*_phasediff.{filetype.value}*",
        f"phasediff {filetype.value} file",
        "",
    )


def get_dwi_fmap_magnitude1_file(filetype: Union[str, DWIFileType]) -> QueryPattern:
    filetype = DWIFileType(filetype)
    return QueryPattern(
        f"fmap/sub-*_ses-*_magnitude1.{filetype.value}*",
        f"magnitude1 {filetype.value} file",
        "",
    )


def get_t1w_linear(cropped: bool) -> QueryPattern:
    return QueryPattern(
        f"*space-MNI152NLin2009cSym{'_desc-Crop' if cropped else ''}_res-1x1x1_T1w.nii.gz",
        (
            "T1w image registered in MNI152NLin2009cSym space "
            f"{'and cropped (matrix size 169×208×179) ' if cropped else ''} "
            "using t1-linear pipeline"
        ),
        "t1-linear",
    )


def get_t1_freesurfer_white_matter_surface(
    hemisphere: Union[str, HemiSphere],
) -> QueryPattern:
    """Return the pattern to query white matter border surface files from the Freesurfer output.

    Parameters
    ----------
    hemisphere : str or HemiSphere
        The hemisphere for which to get the surface.

    Returns
    -------
    Query :
        The query to use with a file reader.
    """
    hemisphere = HemiSphere(hemisphere)
    return QueryPattern(
        f"t1/freesurfer_cross_sectional/sub-*_ses-*/surf/{hemisphere.value}.white",
        (
            f"{'right' if hemisphere == HemiSphere.RIGHT else 'left'} white matter/gray "
            f"matter border surface ({hemisphere.value}.white)."
        ),
        "t1-freesurfer",
    )


def get_t1_freesurfer_longitudinal_white_matter_surface(
    hemisphere: Union[str, HemiSphere],
) -> QueryPattern:
    """Return the query to get white matter border surface files from the Freesurfer longitudinal output.

    Parameters
    ----------
    hemisphere : str or HemiSphere
        The hemisphere for which to get the surface.

    Returns
    -------
    QueryPattern :
        The pattern to use with a file reader.
    """
    hemisphere = HemiSphere(hemisphere)
    return QueryPattern(
        f"t1/long-*/freesurfer_longitudinal/sub-*_ses-*.long.sub-*_*/surf/{hemisphere.value}.white",
        (
            f"{'right' if hemisphere == HemiSphere.RIGHT else 'left'} white matter/gray matter border "
            f"surface ({hemisphere.value}.white) generated with t1-freesurfer-longitudinal."
        ),
        "t1-freesurfer and t1-freesurfer longitudinal",
    )


def _get_annot_file_name(hemisphere: HemiSphere, parcellation: Parcellation) -> str:
    if parcellation == Parcellation.DESIKAN:
        return f"{hemisphere.value}.aparc.annot"
    if parcellation == Parcellation.DESTRIEUX:
        return f"{hemisphere.value}.aparc.a2009s.annot"


def get_t1_freesurfer_segmentation(parcellation: Parcellation) -> QueryPattern:
    parcellation = Parcellation(parcellation)
    filename = (
        f"aparc{'.a2009s' if parcellation == Parcellation.DESTRIEUX else ''}+aseg.mgz"
    )
    return QueryPattern(
        f"t1/freesurfer_cross_sectional/sub-*_ses-*/mri/{filename}",
        f"{parcellation.value}-based segmentation (mri/{filename}).",
        "t1-freesurfer",
    )


def get_t1_freesurfer_parcellation(
    hemisphere: Union[str, HemiSphere],
    parcellation: Union[str, Parcellation],
) -> QueryPattern:
    hemisphere = HemiSphere(hemisphere)
    parcellation = Parcellation(parcellation)
    return QueryPattern(
        f"t1/freesurfer_cross_sectional/sub-*_ses-*/label/{_get_annot_file_name(hemisphere, parcellation)}",
        (
            f"{'left' if hemisphere == HemiSphere.LEFT else 'right'} hemisphere surface-based "
            f"{parcellation.value} parcellation (label/{_get_annot_file_name(hemisphere, parcellation)})."
        ),
        "t1-freesurfer",
    )


def get_t1_freesurfer_template(parcellation: Union[str, Parcellation]) -> QueryPattern:
    parcellation = Parcellation(parcellation)
    filename = (
        f"aparc{'.a2009s' if parcellation == Parcellation.DESTRIEUX else ''}+aseg.mgz"
    )
    return QueryPattern(
        f"freesurfer_unbiased_template/sub-*_long-*/mri/{filename}",
        f"{parcellation.value}-based segmentation (mri/{filename}) from unbiased template.",
        "t1-freesurfer-longitudinal or t1-freesurfer-template",
    )


def get_t1_freesurfer_longitudinal_parcellation(
    hemisphere: Union[str, HemiSphere],
    parcellation: Union[str, Parcellation],
) -> QueryPattern:
    hemisphere = HemiSphere(hemisphere)
    parcellation = Parcellation(parcellation)
    return QueryPattern(
        f"t1/long-*/freesurfer_longitudinal/sub-*_ses-*.long.sub-*_*/label/{_get_annot_file_name(hemisphere, parcellation)}",
        (
            f"{'left' if hemisphere == HemiSphere.LEFT else 'right'} hemisphere surface-based "
            f"{parcellation.value} parcellation (label/{_get_annot_file_name(hemisphere, parcellation)}) "
            "generated with t1-freesurfer-longitudinal."
        ),
        "t1-freesurfer and t1-freesurfer-longitudinal",
    )


def aggregator(func):
    """If the decorated function receives iterable arguments,
    this decorator will call the decorated function for each
    value in the iterable and aggregate the results in a list.
    This works only if the iterables provided have the same length.
    Arguments lefts as non-iterable will be repeated.

    Examples
    --------
    The function `t1_volume_native_tpm` expects an integer defining the
    mask tissue and returns a dictionary describing the files to read:

    >>> import json
    >>> from clinica.utils.input_files import t1_volume_native_tpm
    >>> print(json.dumps(t1_volume_native_tpm(1), indent=3))
    {
        "pattern": "t1/spm/segmentation/native_space/*_*_T1w_segm-graymatter_probability.nii*",
        "description": "Tissue probability map graymatter in native space",
        "needed_pipeline": "t1-volume-tissue-segmentation"
    }

    Without the `aggregator` decorator, querying files for multiple
    tissues would have to be implemented in a loop:

    >>> print(json.dumps([t1_volume_native_tpm(tissue) for tissue in (1, 2)], indent=3))
    [
        {
            "pattern": "t1/spm/segmentation/native_space/*_*_T1w_segm-graymatter_probability.nii*",
            "description": "Tissue probability map graymatter in native space",
            "needed_pipeline": "t1-volume-tissue-segmentation"
        },
        {
            "pattern": "t1/spm/segmentation/native_space/*_*_T1w_segm-whitematter_probability.nii*",
            "description": "Tissue probability map whitematter in native space",
            "needed_pipeline": "t1-volume-tissue-segmentation"
        }
    ]

    Although this is fine, you might not know in a pipeline what was provided (scalar or iterable).
    With the `aggregator` decorator, you can pass both:

    >>> t1_volume_native_tpm((1, 2))
    [
        {
            "pattern": "t1/spm/segmentation/native_space/*_*_T1w_segm-graymatter_probability.nii*",
            "description": "Tissue probability map graymatter in native space",
            "needed_pipeline": "t1-volume-tissue-segmentation"
        },
        {
            "pattern": "t1/spm/segmentation/native_space/*_*_T1w_segm-whitematter_probability.nii*",
            "description": "Tissue probability map whitematter in native space",
            "needed_pipeline": "t1-volume-tissue-segmentation"
        }
    ]

    This works also with multiple args and kwargs:

    >>> from clinica.utils.input_files import t1_volume_native_tpm_in_mni
    >>> print(json.dumps(t1_volume_native_tpm_in_mni(1, False), indent=3))
    {
        "pattern": "t1/spm/segmentation/normalized_space/*_*_T1w_segm-graymatter_space-Ixi549Space_modulated-off_probability.nii*",
        "description": "Tissue probability map graymatter based on native MRI in MNI space (Ixi549) without modulation.",
        "needed_pipeline": "t1-volume-tissue-segmentation"
    }
    >>> print(json.dumps(t1_volume_native_tpm_in_mni(1, (True, False)), indent=3))
    [
        {
            "pattern": "t1/spm/segmentation/normalized_space/*_*_T1w_segm-graymatter_space-Ixi549Space_modulated-on_probability.nii*",
            "description": "Tissue probability map graymatter based on native MRI in MNI space (Ixi549) with modulation.",
            "needed_pipeline": "t1-volume-tissue-segmentation"
        },
        {
            "pattern": "t1/spm/segmentation/normalized_space/*_*_T1w_segm-graymatter_space-Ixi549Space_modulated-off_probability.nii*",
            "description": "Tissue probability map graymatter based on native MRI in MNI space (Ixi549) without modulation.",
            "needed_pipeline": "t1-volume-tissue-segmentation"
        }
    ]
    >>> print(json.dumps(t1_volume_native_tpm_in_mni((1, 2), (True, False)), indent=3))
    [
        {
            "pattern": "t1/spm/segmentation/normalized_space/*_*_T1w_segm-graymatter_space-Ixi549Space_modulated-on_probability.nii*",
            "description": "Tissue probability map graymatter based on native MRI in MNI space (Ixi549) with modulation.",
            "needed_pipeline": "t1-volume-tissue-segmentation"
        },
        {
            "pattern": "t1/spm/segmentation/normalized_space/*_*_T1w_segm-whitematter_space-Ixi549Space_modulated-off_probability.nii*",
            "description": "Tissue probability map whitematter based on native MRI in MNI space (Ixi549) without modulation.",
            "needed_pipeline": "t1-volume-tissue-segmentation"
        }
    ]
    """

    @functools.wraps(func)
    def wrapper_aggregator(*args, **kwargs):
        # Get the lengths of iterable args and kwargs
        arg_sizes = [
            len(arg)
            for arg in args
            if (isinstance(arg, Iterable) and not isinstance(arg, str))
        ]
        arg_sizes += [
            len(arg)
            for k, arg in kwargs.items()
            if (isinstance(arg, Iterable) and not isinstance(arg, str))
        ]

        # If iterable args/kwargs have different lengths, raise
        if len(set(arg_sizes)) > 1:
            raise ValueError(f"Arguments must have the same length.")

        # No iterable case, just call the function
        if len(arg_sizes) == 0:
            return func(*args, **kwargs)

        # Handle args first by repeating non-iterable values
        arg_size = arg_sizes[0]
        new_args = []
        for arg in args:
            if not (isinstance(arg, Iterable) and not isinstance(arg, str)):
                new_args.append((arg,) * arg_size)
            else:
                new_args.append(arg)

        # Same thing for kwargs
        new_kwargs = [{} for _ in range(arg_size)]
        for k, arg in kwargs.items():
            for i in range(len(new_kwargs)):
                if not (isinstance(arg, Iterable) and not isinstance(arg, str)):
                    new_kwargs[i][k] = arg
                else:
                    new_kwargs[i][k] = arg[i]

        # Properly encapsulate in a for loop
        if len(new_args) == 0:
            return [func(**x) for x in new_kwargs]
        elif len(new_kwargs) == 0:
            return [func(*x) for x in zip(*new_args)]
        return [func(*x, **y) for x, y in zip(zip(*new_args), new_kwargs)]

    return wrapper_aggregator


@aggregator
def get_t1_volume_native_tpm(
    tissue_number: int, modulation: bool, mni_space: bool
) -> QueryPattern:
    tissue = get_spm_tissue_from_index(tissue_number)
    description = f"Tissue probability map {tissue.value} "
    pattern_modulation = ""
    space = ""
    if mni_space:
        pattern_modulation = f"_modulated-{'on' if modulation else 'off'}"
        space = "_space-Ixi549Space"
        description += f"based on native MRI in MNI space (Ixi549) {'with' if modulation else 'without'} modulation."
    else:
        description += "in native space"

    return QueryPattern(
        str(
            Path("t1")
            / "spm"
            / "segmentation"
            / f"{'normalized' if mni_space else 'native'}_space"
            / f"*_*_T1w_segm-{tissue.value}{space}{pattern_modulation}_probability.nii*"
        ),
        description,
        "t1-volume-tissue-segmentation",
    )


@aggregator
def get_t1_volume_dartel_input_tissue(tissue_number: int) -> QueryPattern:
    tissue = get_spm_tissue_from_index(tissue_number)
    return QueryPattern(
        str(
            Path("t1")
            / "spm"
            / "segmentation"
            / "dartel_input"
            / f"*_*_T1w_segm-{tissue.value}_dartelinput.nii*"
        ),
        f"Dartel input for tissue probability map {tissue.value} from T1w MRI",
        "t1-volume-tissue-segmentation",
    )


def get_t1_volume_template_tpm_in_mni(
    group_label: str, tissue_number: int, modulation: bool, fwhm: Optional[int] = None
) -> QueryPattern:
    """Build the pattern required by clinica_file_reader to get the tissue
    probability maps based on group template in MNI space.

    Parameters
    ----------
    group_label : str
        Label used for the group of interest.

    tissue_number : int
        An integer defining the tissue of interest.

    modulation : {"on", "off"}
        Whether modulation is on or off.

    fwhm : int, optional
        The smoothing kernel in millimeters.

    Returns
    -------
    QueryPattern :
        Pattern to be passed to clinica_file_reader.
    """
    from .spm import get_spm_tissue_from_index

    tissue = get_spm_tissue_from_index(tissue_number)
    pattern_modulation = "on" if modulation else "off"
    description_modulation = "with" if modulation else "without"
    fwhm_key_value = f"_fwhm-{fwhm}mm" if fwhm else ""
    fwhm_description = f"with {fwhm}mm smoothing" if fwhm else "with no smoothing"

    return QueryPattern(
        str(
            Path("t1")
            / "spm"
            / "dartel"
            / f"group-{group_label}"
            / f"*_T1w_segm-{tissue.value}_space-Ixi549Space_modulated-{pattern_modulation}{fwhm_key_value}_probability.nii*"
        ),
        (
            f"Tissue probability map {tissue.value} based on {group_label} template in MNI space "
            f"(Ixi549) {description_modulation} modulation and {fwhm_description}."
        ),
        "t1-volume",
    )


def get_t1_volume_deformation_to_template(group_label: str) -> QueryPattern:
    return QueryPattern(
        str(
            Path("t1")
            / "spm"
            / "dartel"
            / f"group-{group_label}"
            / f"sub-*_ses-*_T1w_target-{group_label}_transformation-forward_deformation.nii*"
        ),
        f"Deformation from native space to group template {group_label} space.",
        "t1-volume-create-dartel",
    )


@aggregator
def get_t1_volume_i_th_iteration_group_template(
    group_label: str, i: int
) -> QueryPattern:
    return QueryPattern(
        str(
            Path(f"group-{group_label}")
            / "t1"
            / f"group-{group_label}_iteration-{i}_template.nii*"
        ),
        f"Iteration #{i} of Dartel template {group_label}",
        "t1-volume or t1-volume-create-dartel",
    )


def get_t1_volume_group_template(group_label: str) -> QueryPattern:
    return QueryPattern(
        str(Path(f"group-{group_label}") / "t1" / f"group-{group_label}_template.nii*"),
        f"T1w template file of group {group_label}",
        "t1-volume or t1-volume-create-dartel",
    )


def get_dwi_dti(
    measure: Union[str, DTIBasedMeasure], space: Optional[str] = None
) -> QueryPattern:
    """Return the query pattern required to capture DWI DTI images.

    Parameters
    ----------
    measure : DTIBasedMeasure or str
        The DTI based measure to consider.

    space : str, optional
        The space to consider.
        By default, all spaces are considered (i.e. '*' is used in regexp).

    Returns
    -------
    QueryPattern :
        The query pattern to get DWI DTI images.
    """
    measure = DTIBasedMeasure(measure)
    space = space or "*"

    return QueryPattern(
        f"dwi/dti_based_processing/*/*_space-{space}_{measure.value}.nii.gz",
        f"DTI-based {measure.value} in space {space}.",
        "dwi_dti",
    )


def get_pet_nifti(
    tracer: Optional[Union[str, Tracer]] = None,
    reconstruction: Optional[Union[str, ReconstructionMethod]] = None,
) -> QueryPattern:
    """Return the query pattern required to capture PET scans.

    Parameters
    ----------
    tracer : Tracer, optional
        If specified, the query will only match PET scans acquired
        with the requested tracer.
        If None, the query will match all PET sans independently of
        the tracer used.

    reconstruction : ReconstructionMethod, optional
        If specified, the query will only match PET scans reconstructed
        with the specified method.
        If None, the query will match all PET scans independently of the
        reconstruction method used.

    Returns
    -------
    QueryPattern :
        The query pattern to get PET scans.
    """
    description = f"PET data"
    trc = ""
    if tracer is not None:
        tracer = Tracer(tracer)
        trc = f"_trc-{tracer.value}"
        description += f" with {tracer.value} tracer"
    rec = ""
    if reconstruction is not None:
        reconstruction = ReconstructionMethod(reconstruction)
        rec = f"_rec-{reconstruction.value}"
        description += f" and reconstruction method {reconstruction.value}"

    return QueryPattern(
        str(Path("pet") / f"*{trc}{rec}_pet.nii*"),
        description,
        "",
    )


def get_pet_volume_normalized_suvr(
    tracer: Union[str, Tracer],
    group_label: str,
    suvr_reference_region: Union[str, SUVRReferenceRegion],
    use_brainmasked_image: bool,
    use_pvc_data: bool,
    fwhm: int = 0,
) -> QueryPattern:
    tracer = Tracer(tracer)
    region = SUVRReferenceRegion(suvr_reference_region)

    if use_brainmasked_image:
        mask_key_value = "_mask-brain"
        mask_description = "brain-masked"
    else:
        mask_key_value = ""
        mask_description = "full"
    if use_pvc_data:
        pvc_key_value = "_pvc-rbv"
        pvc_description = "using RBV method for PVC"
    else:
        pvc_key_value = ""
        pvc_description = "without PVC"
    if fwhm:
        fwhm_key_value = f"_fwhm-{fwhm}mm"
        fwhm_description = f"with {fwhm}mm smoothing"
    else:
        fwhm_key_value = ""
        fwhm_description = "with no smoothing"
    suvr_key_value = f"_suvr-{region.value}"

    return QueryPattern(
        str(
            Path("pet")
            / "preprocessing"
            / f"group-{group_label}"
            / f"*_trc-{tracer.value}_pet_space-Ixi549Space{pvc_key_value}{suvr_key_value}{mask_key_value}{fwhm_key_value}_pet.nii*"
        ),
        (
            f"{mask_description} SUVR map (using {region.value} region) of {tracer.value}-PET "
            f"{pvc_description} and {fwhm_description} in Ixi549Space space based on {group_label} DARTEL template"
        ),
        "pet-volume",
    )


def get_pet_linear_nifti(
    tracer: Union[str, Tracer],
    suvr_reference_region: Union[str, SUVRReferenceRegion],
    uncropped_image: bool,
) -> QueryPattern:
    tracer = Tracer(tracer)
    region = SUVRReferenceRegion(suvr_reference_region)
    description = "" if uncropped_image else "_desc-Crop"

    return QueryPattern(
        str(
            Path("pet_linear")
            / f"*_trc-{tracer.value}_pet_space-MNI152NLin2009cSym{description}_res-1x1x1_suvr-{region.value}_pet.nii.gz"
        ),
        "PET nifti image obtained with pet-linear",
        "pet-linear",
    )
