"""This module contains dictionaries used in inputs.py::clinica_{file|group}_reader().

These dictionaries describe files to grab.
"""

import functools
from collections import UserString
from collections.abc import Iterable
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Optional, Union

from clinica.utils.dwi import DTIBasedMeasure
from clinica.utils.group import GroupID
from clinica.utils.image import HemiSphere
from clinica.utils.pet import ReconstructionMethod, SUVRReferenceRegion, Tracer

from .spm import SPMTissue, get_spm_tissue_from_index

__all__ = [
    "DWIFileType",
    "Parcellation",
    "QueryPattern",
    "get_t1w_mri",
    "get_t2w_mri",
    "get_t1_freesurfer_segmentation_white_matter",
    "get_t1_freesurfer_extracted_brain",
    "get_t1_freesurfer_intensity_normalized_volume_after_nu",
    "get_t1_freesurfer_longitudinal_intensity_normalized_volume_after_nu",
    "get_t1w_to_mni_transform",
    "get_dwi_file",
    "get_dwi_preprocessed_file",
    "get_dwi_preprocessed_brainmask",
    "get_dwi_fmap_phasediff_file",
    "get_dwi_fmap_magnitude1_file",
    "get_t1w_linear",
    "get_t1_freesurfer_white_matter_surface",
    "get_t1_freesurfer_longitudinal_white_matter_surface",
    "get_t1_freesurfer_segmentation",
    "get_t1_freesurfer_statistics",
    "get_t1_freesurfer_parcellation",
    "get_t1_freesurfer_template",
    "get_t1_freesurfer_longitudinal_parcellation",
    "get_t1_volume_tpm",
    "get_t1_volume_dartel_input_tissue",
    "get_t1_volume_template_tpm_in_mni",
    "get_t1_volume_deformation_to_template",
    "get_t1_volume_i_th_iteration_group_template",
    "get_t1_volume_group_template",
    "get_dwi_dti",
    "get_pet_nifti",
    "get_pet_volume_normalized_suvr",
    "get_pet_linear_nifti",
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


class Parcellation(str, Enum):
    """The possible atlas names used for deriving parcellations and segmentations."""

    DESIKAN = "Desikan"
    DESTRIEUX = "Destrieux"


class DWIFileType(str, Enum):
    NII = "nii"
    JSON = "json"
    BVEC = "bvec"
    BVAL = "bval"


def get_t1w_mri() -> QueryPattern:
    """Get T1W MRI in BIDS."""
    return QueryPattern("sub-*_ses-*_t1w.nii*", "T1w MRI", "")


def get_t2w_mri() -> QueryPattern:
    """Get T2W FLAIR MRI in BIDS."""
    return QueryPattern("sub-*_ses-*_flair.nii*", "FLAIR T2w MRI", "")


def get_t1_freesurfer_segmentation_white_matter() -> QueryPattern:
    """GET Freesurfer segmentation of white matter."""
    return QueryPattern(
        "t1/freesurfer_cross_sectional/sub-*_ses-*/mri/wm.seg.mgz",
        "segmentation of white matter (mri/wm.seg.mgz).",
        "t1-freesurfer",
    )


def get_t1_freesurfer_extracted_brain() -> QueryPattern:
    return QueryPattern(
        "t1/freesurfer_cross_sectional/sub-*_ses-*/mri/brain.mgz",
        "extracted brain from T1w MRI (mri/brain.mgz).",
        "t1-freesurfer",
    )


def get_t1_freesurfer_intensity_normalized_volume_after_nu() -> QueryPattern:
    return QueryPattern(
        "t1/freesurfer_cross_sectional/sub-*_ses-*/mri/orig_nu.mgz",
        (
            "intensity normalized volume generated after correction for "
            "non-uniformity in FreeSurfer (mri/orig_nu.mgz)."
        ),
        "t1-freesurfer",
    )


def get_t1_freesurfer_longitudinal_intensity_normalized_volume_after_nu() -> (
    QueryPattern
):
    return QueryPattern(
        "t1/long-*/freesurfer_longitudinal/sub-*_ses-*.long.sub-*_*/mri/orig_nu.mgz",
        (
            "intensity normalized volume generated after correction for "
            "non-uniformity in FreeSurfer (orig_nu.mgz) in longitudinal"
        ),
        "t1-freesurfer and t1-freesurfer longitudinal",
    )


def get_t1w_to_mni_transform() -> QueryPattern:
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


def get_dwi_preprocessed_brainmask() -> QueryPattern:
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


def get_t1_freesurfer_statistics(
    atlas: str, hemisphere: Union[str, HemiSphere]
) -> QueryPattern:
    hemisphere = HemiSphere(hemisphere)
    return QueryPattern(
        f"t1/freesurfer_cross_sectional/sub-*_ses-*/stats/{hemisphere.value}.{atlas}.stats",
        f"{atlas}-based segmentation",
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
            if (isinstance(arg, Iterable) and not isinstance(arg, (str, UserString)))
        ]
        arg_sizes += [
            len(arg)
            for k, arg in kwargs.items()
            if (isinstance(arg, Iterable) and not isinstance(arg, (str, UserString)))
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
            if not (
                isinstance(arg, Iterable) and not isinstance(arg, (str, UserString))
            ):
                new_args.append((arg,) * arg_size)
            else:
                new_args.append(arg)

        # Same thing for kwargs
        new_kwargs = [{} for _ in range(arg_size)]
        for k, arg in kwargs.items():
            for i in range(len(new_kwargs)):
                if not (
                    isinstance(arg, Iterable) and not isinstance(arg, (str, UserString))
                ):
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
def get_t1_volume_tpm(
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
    group_id: GroupID,
    tissue: Union[int, SPMTissue],
    modulation: bool,
    fwhm: Optional[int] = None,
) -> QueryPattern:
    """Build the pattern required by clinica_file_reader to get the tissue
    probability maps based on group template in MNI space.

    Parameters
    ----------
    group_id : GroupID
        The ID for the group of interest.

    tissue : int or SPMTissue
        Either the tissue of interest, or an integer defining the tissue of interest.

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

    tissue = get_spm_tissue_from_index(tissue)
    pattern_modulation = "on" if modulation else "off"
    description_modulation = "with" if modulation else "without"
    fwhm_key_value = ""
    fwhm_description = "with no smoothing"
    if fwhm is not None and fwhm != 0:
        fwhm_key_value = f"_fwhm-{fwhm}mm"
        fwhm_description = f"with {fwhm}mm smoothing"

    return QueryPattern(
        str(
            Path("t1")
            / "spm"
            / "dartel"
            / str(group_id)
            / f"*_T1w_segm-{tissue.value}_space-Ixi549Space_modulated-{pattern_modulation}{fwhm_key_value}_probability.nii*"
        ),
        (
            f"Tissue probability map {tissue.value} based on {group_id.label} template in MNI space "
            f"(Ixi549) {description_modulation} modulation and {fwhm_description}."
        ),
        "t1-volume",
    )


def get_t1_volume_deformation_to_template(group_id: GroupID) -> QueryPattern:
    return QueryPattern(
        str(
            Path("t1")
            / "spm"
            / "dartel"
            / str(group_id)
            / f"sub-*_ses-*_T1w_target-{group_id.label}_transformation-forward_deformation.nii*"
        ),
        f"Deformation from native space to group template {group_id.label} space.",
        "t1-volume-create-dartel",
    )


@aggregator
def get_t1_volume_i_th_iteration_group_template(
    group_id: GroupID, i: int
) -> QueryPattern:
    return QueryPattern(
        str(Path(str(group_id)) / "t1" / f"{group_id}_iteration-{i}_template.nii*"),
        f"Iteration #{i} of Dartel template {group_id.label}",
        "t1-volume or t1-volume-create-dartel",
    )


def get_t1_volume_group_template(group_id: GroupID) -> QueryPattern:
    return QueryPattern(
        str(Path(str(group_id)) / "t1" / f"{group_id}_template.nii*"),
        f"T1w template file of group {group_id.label}",
        "t1-volume or t1-volume-create-dartel",
    )


def get_dwi_dti(
    measure: Union[str, DTIBasedMeasure],
    space: Optional[str] = None,
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
    tracer: Optional[Union[str, Tracer]] = None,
    group_id: Optional[GroupID] = None,
    suvr_reference_region: Optional[Union[str, SUVRReferenceRegion]] = None,
    use_brainmasked_image: bool = False,
    use_pvc_data: bool = False,
    fwhm: int = 0,
) -> QueryPattern:
    group_description = "any available template found"
    if group_id:
        group_description = f"{group_id.label} DARTEL template"
    else:
        group_id = "*"
    tracer_key_value = ""
    tracer_description = ""
    if tracer:
        tracer = Tracer(tracer)
        tracer_key_value = f"_trc-{tracer.value}"
        tracer_description = f" of {tracer.value}-PET"
    suvr_key_value = ""
    region_description = ""
    if suvr_reference_region:
        region = SUVRReferenceRegion(suvr_reference_region)
        suvr_key_value = f"_suvr-{region.value}"
        region_description = f" (using {region.value} region)"
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

    return QueryPattern(
        str(
            Path("pet")
            / "preprocessing"
            / str(group_id)
            / f"*{tracer_key_value}_pet_space-Ixi549Space{pvc_key_value}{suvr_key_value}{mask_key_value}{fwhm_key_value}_pet.nii*"
        ),
        (
            f"{mask_description} SUVR map{region_description}{tracer_description} "
            f"{pvc_description} and {fwhm_description} in Ixi549Space space based on {group_description}"
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
