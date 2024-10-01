"""This module contains dictionaries used in inputs.py::clinica_{file|group}_reader().

These dictionaries describe files to grab.
"""

import functools
from collections.abc import Iterable
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Optional, Union

from clinica.utils.dwi import DTIBasedMeasure
from clinica.utils.image import HemiSphere
from clinica.utils.pet import ReconstructionMethod, SUVRReferenceRegion, Tracer

__all__ = [
    "Query",
    "QueryName",
    "query_factory",
]


@dataclass
class Query:
    """Represents a query for the clinica_file_reader.

    Attributes
    ----------
    pattern : str
        The pattern used to match file names.

    description : str
        A plain text description of the files the query matches.

    needed_pipeline : list of tuples of str
        The pipelines that should have been run in order to have the requested files.
    """

    pattern: str
    description: str
    needed_pipeline: str


class QueryName(str, Enum):
    """The different names for usual queries in Clinica.

    T1W : Get T1W MRI in BIDS
    T2W : Get T2W FLAIR MRI in BIDS
    T1_FS_WM : GET Freesurfer segmentation of white matter
    T1_FS_BRAIN :  Get Freesurfer extracted brain from T1w MRI
    T1_FS_ORIG_NU : Get Freesurfer intensity normalized volume after correction for non-uniformity
    T1_FS_WM_SURF : Get white matter border surface files from the Freesurfer output
    T1_FS_LONG_SURF : Get white matter border surface files from the Freesurfer longitudinal output
    """

    T1W = "T1W"
    T2W = "T2W"
    T1_FS_WM = "T1_FS_WM"
    T1_FS_BRAIN = "T1_FS_BRAIN"
    T1_FS_ORIG_NU = "T1_FS_ORIG_NU"
    T1_FS_LONG_ORIG_NU = "T1_FS_LONG_ORIG_NU"
    T1_FS_WM_SURF = "T1_FS_WM_SURF"
    T1_FS_LONG_SURF = "T1_FS_LONG_SURF"
    T1W_LINEAR = "T1W_LINEAR"
    T1W_TO_MNI_TRANSFORM = "T1W_TO_MNI_TRANSFORM"
    T1_FS_PARC = "T1_FS_PARC"
    T1_FS_LONG_PARC = "T1_FS_LONG_PARC"
    T1_FS_SEG = "T1_FS_SEG"
    T1_FS_TEMPLATE = "T1_FS_TEMPLATE"
    DWI = "DWI"
    DWI_PREPROC = "DWI_PREPROC"
    DWI_PREPROC_BRAINMASK = "DWI_PREPROC_BRAINMASK"
    DWI_FMAP_PHASEDIFF = "DWI_FMAP_PHASEDIFF"
    DWI_FMAP_MAGNITUDE1 = "DWI_FMAP_MAGNITUDE1"
    DWI_DTI = "DWI_DTI"


class Parcellation(str, Enum):
    DESIKAN = "Desikan"
    DESTRIEUX = "Destrieux"


class DWIFileType(str, Enum):
    NII = "nii"
    JSON = "json"
    BVEC = "bvec"
    BVAL = "bval"


def query_factory(name: Union[str, QueryName], *args, **kwargs) -> Query:
    """Return the query corresponding to the provided name.

    Additional arguments can be passed if the query builder is parametric.

    Parameters
    ----------
    name : str or QueryName
        The name of the desired query.

    Returns
    -------
    Query :
        The desired query.
    """
    name = QueryName(name)
    if name == QueryName.T1W:
        return Query("sub-*_ses-*_t1w.nii*", "T1w MRI", "")
    if name == QueryName.T2W:
        return Query("sub-*_ses-*_flair.nii*", "FLAIR T2w MRI", "")
    if name == QueryName.T1_FS_WM:
        return Query(
            "t1/freesurfer_cross_sectional/sub-*_ses-*/mri/wm.seg.mgz",
            "segmentation of white matter (mri/wm.seg.mgz).",
            "t1-freesurfer",
        )
    if name == QueryName.T1_FS_BRAIN:
        return Query(
            "t1/freesurfer_cross_sectional/sub-*_ses-*/mri/brain.mgz",
            "extracted brain from T1w MRI (mri/brain.mgz).",
            "t1-freesurfer",
        )
    if name == QueryName.T1_FS_ORIG_NU:
        return Query(
            "t1/freesurfer_cross_sectional/sub-*_ses-*/mri/orig_nu.mgz",
            (
                "intensity normalized volume generated after correction for "
                "non-uniformity in FreeSurfer (mri/orig_nu.mgz)."
            ),
            "t1-freesurfer",
        )
    if name == QueryName.T1_FS_LONG_ORIG_NU:
        return Query(
            "t1/long-*/freesurfer_longitudinal/sub-*_ses-*.long.sub-*_*/mri/orig_nu.mgz",
            (
                "intensity normalized volume generated after correction for "
                "non-uniformity in FreeSurfer (orig_nu.mgz) in longitudinal"
            ),
            "t1-freesurfer and t1-freesurfer longitudinal",
        )
    if name == QueryName.T1_FS_WM:
        return t1_freesurfer_white_matter_surface(*args, **kwargs)
    if name == QueryName.T1_FS_LONG_SURF:
        return t1_freesurfer_longitudinal_white_matter_surface(*args, **kwargs)
    if name == QueryName.T1W_LINEAR:
        return get_t1w_linear(*args, **kwargs)
    if name == QueryName.T1W_TO_MNI_TRANSFORM:
        return Query(
            "*space-MNI152NLin2009cSym_res-1x1x1_affine.mat",
            "Transformation matrix from T1W image to MNI space using t1-linear pipeline",
            "t1-linear",
        )
    if name == QueryName.T1_FS_PARC:
        return get_t1_freesurfer_parcellation(*args, **kwargs)
    if name == QueryName.T1_FS_LONG_PARC:
        return get_t1_freesurfer_longitudinal_parcellation(*args, **kwargs)
    if name == QueryName.T1_FS_SEG:
        return get_t1_freesurfer_segmentation(*args, **kwargs)
    if name == QueryName.T1_FS_TEMPLATE:
        return get_t1_freesurfer_template(*args, **kwargs)
    if name == QueryName.DWI:
        return get_dwi_file(*args, **kwargs)
    if name == QueryName.DWI_PREPROC:
        return get_dwi_preprocessed_file(*args, **kwargs)
    if name == QueryName.DWI_PREPROC_BRAINMASK:
        return Query(
            "dwi/preprocessing/sub-*_ses-*_space-*_brainmask.nii*",
            "b0 brainmask",
            "dwi-preprocessing-using-t1 or dwi-preprocessing-using-fieldmap",
        )
    if name == QueryName.DWI_FMAP_PHASEDIFF:
        return get_dwi_fmap_phasediff_file(*args, **kwargs)
    if name == QueryName.DWI_FMAP_MAGNITUDE1:
        return get_dwi_fmap_magnitude1_file(*args, **kwargs)
    if name == QueryName.DWI_DTI:
        return dwi_dti(*args, **kwargs)


def get_dwi_file(filetype: Union[str, DWIFileType]) -> Query:
    filetype = DWIFileType(filetype)
    return Query(
        f"dwi/sub-*_ses-*_dwi.{filetype.value}*", f"DWI {filetype.value} files.", ""
    )


def get_dwi_preprocessed_file(filetype: Union[str, DWIFileType]) -> Query:
    filetype = DWIFileType(filetype)
    return Query(
        f"dwi/preprocessing/sub-*_ses-*_space-*_desc-preproc_dwi.{filetype.value}*",
        f"preprocessed {filetype.value} files",
        "dwi-preprocessing-using-t1 or dwi-preprocessing-using-fieldmap",
    )


def get_dwi_fmap_phasediff_file(filetype: Union[str, DWIFileType]) -> Query:
    filetype = DWIFileType(filetype)
    return Query(
        f"fmap/sub-*_ses-*_phasediff.{filetype.value}",
        f"phasediff {filetype.value} file",
        "",
    )


def get_dwi_fmap_magnitude1_file(filetype: Union[str, DWIFileType]) -> Query:
    filetype = DWIFileType(filetype)
    return Query(
        f"fmap/sub-*_ses-*_magnitude1.{filetype.value}*",
        f"magnitude1 {filetype.value} file",
        "",
    )


def get_t1w_linear(cropped: bool) -> Query:
    return Query(
        f"*space-MNI152NLin2009cSym{'_desc-Crop' if cropped else ''}_res-1x1x1_T1w.nii.gz",
        (
            "T1w image registered in MNI152NLin2009cSym space "
            f"{'and cropped (matrix size 169×208×179) ' if cropped else ''} "
            "using t1-linear pipeline"
        ),
        "t1-linear",
    )


def t1_freesurfer_white_matter_surface(hemisphere: Union[str, HemiSphere]) -> Query:
    """Return the query to get white matter border surface files from the Freesurfer output.

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
    return Query(
        f"t1/freesurfer_cross_sectional/sub-*_ses-*/surf/{hemisphere.value}.white",
        (
            f"{'right' if hemisphere == HemiSphere.RIGHT else 'left'} white matter/gray "
            f"matter border surface ({hemisphere.value}.white)."
        ),
        "t1-freesurfer",
    )


def t1_freesurfer_longitudinal_white_matter_surface(
    hemisphere: Union[str, HemiSphere],
) -> Query:
    """Return the query to get white matter border surface files from the Freesurfer longitudinal output.

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
    return Query(
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


def get_t1_freesurfer_segmentation(parcellation: Parcellation) -> Query:
    parcellation = Parcellation(parcellation)
    filename = (
        f"aparc{'.a2009s' if parcellation == Parcellation.DESTRIEUX else ''}+aseg.mgz"
    )
    return Query(
        f"t1/freesurfer_cross_sectional/sub-*_ses-*/mri/{filename}",
        f"{parcellation.value}-based segmentation (mri/{filename}).",
        "t1-freesurfer",
    )


def get_t1_freesurfer_parcellation(
    hemisphere: Union[str, HemiSphere],
    parcellation: Union[str, Parcellation],
) -> Query:
    hemisphere = HemiSphere(hemisphere)
    parcellation = Parcellation(parcellation)
    return Query(
        f"t1/freesurfer_cross_sectional/sub-*_ses-*/label/{_get_annot_file_name(hemisphere, parcellation)}",
        (
            f"{'left' if hemisphere == HemiSphere.LEFT else 'right'} hemisphere surface-based "
            f"{parcellation.value} parcellation (label/{_get_annot_file_name(hemisphere, parcellation)})."
        ),
        "t1-freesurfer",
    )


def get_t1_freesurfer_template(parcellation: Parcellation) -> Query:
    parcellation = Parcellation(parcellation)
    filename = (
        f"aparc{'.a2009s' if parcellation == Parcellation.DESTRIEUX else ''}+aseg.mgz"
    )
    return Query(
        f"freesurfer_unbiased_template/sub-*_long-*/mri/{filename}",
        f"{parcellation.value}-based segmentation (mri/{filename}) from unbiased template.",
        "t1-freesurfer-longitudinal or t1-freesurfer-template",
    )


def get_t1_freesurfer_longitudinal_parcellation(
    hemisphere: Union[str, HemiSphere],
    parcellation: Union[str, Parcellation],
) -> Query:
    hemisphere = HemiSphere(hemisphere)
    parcellation = Parcellation(parcellation)
    return Query(
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
def t1_volume_native_tpm(tissue_number: int) -> Query:
    from .spm import get_spm_tissue_from_index

    tissue = get_spm_tissue_from_index(tissue_number)
    return Query(
        str(
            Path("t1")
            / "spm"
            / "segmentation"
            / "native_space"
            / f"*_*_T1w_segm-{tissue.value}_probability.nii*"
        ),
        f"Tissue probability map {tissue.value} in native space",
        "t1-volume-tissue-segmentation",
    )


@aggregator
def t1_volume_dartel_input_tissue(tissue_number: int) -> Query:
    from .spm import get_spm_tissue_from_index

    tissue = get_spm_tissue_from_index(tissue_number)
    return Query(
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


@aggregator
def t1_volume_native_tpm_in_mni(tissue_number: int, modulation: bool) -> Query:
    from .spm import get_spm_tissue_from_index

    tissue = get_spm_tissue_from_index(tissue_number)
    pattern_modulation = "on" if modulation else "off"
    description_modulation = "with" if modulation else "without"

    return Query(
        str(
            Path("t1")
            / "spm"
            / "segmentation"
            / "normalized_space"
            / f"*_*_T1w_segm-{tissue.value}_space-Ixi549Space_modulated-{pattern_modulation}_probability.nii*"
        ),
        (
            f"Tissue probability map {tissue.value} based on "
            f"native MRI in MNI space (Ixi549) {description_modulation} modulation."
        ),
        "t1-volume-tissue-segmentation",
    )


def t1_volume_template_tpm_in_mni(
    group_label: str, tissue_number: int, modulation: bool, fwhm: Optional[int] = None
) -> Query:
    """Build the dictionary required by clinica_file_reader to get the tissue
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
    dict :
        Information dict to be passed to clinica_file_reader.
    """
    from .spm import get_spm_tissue_from_index

    tissue = get_spm_tissue_from_index(tissue_number)
    pattern_modulation = "on" if modulation else "off"
    description_modulation = "with" if modulation else "without"
    fwhm_key_value = f"_fwhm-{fwhm}mm" if fwhm else ""
    fwhm_description = f"with {fwhm}mm smoothing" if fwhm else "with no smoothing"

    return Query(
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


def t1_volume_deformation_to_template(group_label: str) -> Query:
    return Query(
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
def t1_volume_i_th_iteration_group_template(group_label: str, i: int) -> Query:
    return Query(
        str(
            Path(f"group-{group_label}")
            / "t1"
            / f"group-{group_label}_iteration-{i}_template.nii*"
        ),
        f"Iteration #{i} of Dartel template {group_label}",
        "t1-volume or t1-volume-create-dartel",
    )


def t1_volume_final_group_template(group_label: str) -> Query:
    return Query(
        str(Path(f"group-{group_label}") / "t1" / f"group-{group_label}_template.nii*"),
        f"T1w template file of group {group_label}",
        "t1-volume or t1-volume-create-dartel",
    )


def custom_group(pattern, description):
    information = {"pattern": pattern, "description": description}
    return information


def dwi_dti(measure: Union[str, DTIBasedMeasure], space: Optional[str] = None) -> Query:
    """Return the query dict required to capture DWI DTI images.

    Parameters
    ----------
    measure : DTIBasedMeasure or str
        The DTI based measure to consider.

    space : str, optional
        The space to consider.
        By default, all spaces are considered (i.e. '*' is used in regexp).

    Returns
    -------
    dict :
        The query dictionary to get DWI DTI images.
    """
    measure = DTIBasedMeasure(measure)
    space = space or "*"

    return Query(
        f"dwi/dti_based_processing/*/*_space-{space}_{measure.value}.nii.gz",
        f"DTI-based {measure.value} in space {space}.",
        "dwi_dti",
    )


def bids_pet_nii(
    tracer: Optional[Union[str, Tracer]] = None,
    reconstruction: Optional[Union[str, ReconstructionMethod]] = None,
) -> Query:
    """Return the query dict required to capture PET scans.

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
    dict :
        The query dictionary to get PET scans.
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

    return Query(
        str(Path("pet") / f"*{trc}{rec}_pet.nii*"),
        description,
        "",
    )


def pet_volume_normalized_suvr_pet(
    tracer: Union[str, Tracer],
    group_label: str,
    suvr_reference_region: Union[str, SUVRReferenceRegion],
    use_brainmasked_image: bool,
    use_pvc_data: bool,
    fwhm: int = 0,
) -> Query:
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

    return Query(
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


def pet_linear_nii(
    tracer: Union[str, Tracer],
    suvr_reference_region: Union[str, SUVRReferenceRegion],
    uncropped_image: bool,
) -> Query:
    from pathlib import Path

    tracer = Tracer(tracer)
    region = SUVRReferenceRegion(suvr_reference_region)
    description = "" if uncropped_image else "_desc-Crop"

    return Query(
        str(
            Path("pet_linear")
            / f"*_trc-{tracer.value}_pet_space-MNI152NLin2009cSym{description}_res-1x1x1_suvr-{region.value}_pet.nii.gz"
        ),
        "PET nifti image obtained with pet-linear",
        "pet-linear",
    )


# CUSTOM
def custom_pipeline(pattern, description):
    information = {"pattern": pattern, "description": description}
    return information
