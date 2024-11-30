from pathlib import Path

import pytest

from clinica.utils.dwi import DTIBasedMeasure
from clinica.utils.input_files import QueryPattern
from clinica.utils.pet import ReconstructionMethod, Tracer


def test_aggregator():
    from clinica.utils.input_files import aggregator

    @aggregator
    def toy_func(x):
        return x**2

    assert toy_func(2) == 4
    assert toy_func((1, 3, 5)) == [1, 9, 25]

    @aggregator
    def toy_func_2(x, y, z):
        return x + y, z * 3

    assert toy_func_2(1, 2, 3) == (3, 9)
    assert toy_func_2((1, 2), (3, 4), (5, 6)) == [(4, 15), (6, 18)]

    @aggregator
    def toy_func_3(x, y=2, z=3):
        return x * y, x + y + z

    assert toy_func_3(1, 3, 4) == (3, 8)
    assert toy_func_3(1) == (2, 6)
    assert toy_func_3((1, 2)) == [(2, 6), (4, 7)]
    assert toy_func_3((1, 2), y=(3, 4)) == [(3, 7), (8, 9)]
    assert toy_func_3((1, 2), z=(4, 5)) == [(2, 7), (4, 9)]
    assert toy_func_3(1, y=(3, 5)) == [(3, 7), (5, 9)]
    with pytest.raises(
        ValueError,
        match="Arguments must have the same length.",
    ):
        toy_func_3((1, 2, 3), z=(4, 5))


def test_get_t1w_mri():
    from clinica.utils.input_files import get_t1w_mri

    pattern = get_t1w_mri()

    assert pattern.pattern == "sub-*_ses-*_t1w.nii*"
    assert pattern.description == "T1w MRI"
    assert pattern.needed_pipeline == ""


def test_get_t2w_mri():
    from clinica.utils.input_files import get_t2w_mri

    pattern = get_t2w_mri()

    assert pattern.pattern == "sub-*_ses-*_flair.nii*"
    assert pattern.description == "FLAIR T2w MRI"
    assert pattern.needed_pipeline == ""


def test_get_t1_freesurfer_segmentation_white_matter():
    from clinica.utils.input_files import get_t1_freesurfer_segmentation_white_matter

    pattern = get_t1_freesurfer_segmentation_white_matter()

    assert pattern.pattern == "t1/freesurfer_cross_sectional/sub-*_ses-*/mri/wm.seg.mgz"
    assert pattern.description == "segmentation of white matter (mri/wm.seg.mgz)."
    assert pattern.needed_pipeline == "t1-freesurfer"


def test_get_t1_freesurfer_extracted_brain():
    from clinica.utils.input_files import get_t1_freesurfer_extracted_brain

    pattern = get_t1_freesurfer_extracted_brain()

    assert pattern.pattern == "t1/freesurfer_cross_sectional/sub-*_ses-*/mri/brain.mgz"
    assert pattern.description == "extracted brain from T1w MRI (mri/brain.mgz)."
    assert pattern.needed_pipeline == "t1-freesurfer"


def test_get_t1_freesurfer_intensity_normalized_volume_after_nu():
    from clinica.utils.input_files import (
        get_t1_freesurfer_intensity_normalized_volume_after_nu,
    )

    pattern = get_t1_freesurfer_intensity_normalized_volume_after_nu()

    assert (
        pattern.pattern == "t1/freesurfer_cross_sectional/sub-*_ses-*/mri/orig_nu.mgz"
    )
    assert pattern.description == (
        "intensity normalized volume generated after correction for"
        " non-uniformity in FreeSurfer (mri/orig_nu.mgz)."
    )
    assert pattern.needed_pipeline == "t1-freesurfer"


def test_get_t1_freesurfer_longitudinal_intensity_normalized_volume_after_nu():
    from clinica.utils.input_files import (
        get_t1_freesurfer_longitudinal_intensity_normalized_volume_after_nu,
    )

    pattern = get_t1_freesurfer_longitudinal_intensity_normalized_volume_after_nu()

    assert (
        pattern.pattern
        == "t1/long-*/freesurfer_longitudinal/sub-*_ses-*.long.sub-*_*/mri/orig_nu.mgz"
    )
    assert pattern.description == (
        "intensity normalized volume generated after correction for non-uniformity "
        "in FreeSurfer (orig_nu.mgz) in longitudinal"
    )
    assert pattern.needed_pipeline == "t1-freesurfer and t1-freesurfer longitudinal"


def test_get_t1w_to_mni_transform():
    from clinica.utils.input_files import get_t1w_to_mni_transform

    pattern = get_t1w_to_mni_transform()

    assert pattern.pattern == "*space-MNI152NLin2009cSym_res-1x1x1_affine.mat"
    assert (
        pattern.description
        == "Transformation matrix from T1W image to MNI space using t1-linear pipeline"
    )
    assert pattern.needed_pipeline == "t1-linear"


def test_get_dwi_preprocessed_brainmask():
    from clinica.utils.input_files import get_dwi_preprocessed_brainmask

    pattern = get_dwi_preprocessed_brainmask()

    assert pattern.pattern == "dwi/preprocessing/sub-*_ses-*_space-*_brainmask.nii*"
    assert pattern.description == "b0 brainmask"
    assert (
        pattern.needed_pipeline
        == "dwi-preprocessing-using-t1 or dwi-preprocessing-using-fieldmap"
    )


@pytest.mark.parametrize(
    "filetype,expected_pattern,expected_description,expected_pipelines",
    [
        ("nii", "dwi/sub-*_ses-*_dwi.nii*", "DWI nii files.", ""),
        ("json", "dwi/sub-*_ses-*_dwi.json*", "DWI json files.", ""),
        ("bvec", "dwi/sub-*_ses-*_dwi.bvec*", "DWI bvec files.", ""),
        ("bval", "dwi/sub-*_ses-*_dwi.bval*", "DWI bval files.", ""),
    ],
)
def test_get_dwi_file(
    filetype: str,
    expected_pattern: str,
    expected_description: str,
    expected_pipelines: str,
):
    from clinica.utils.input_files import get_dwi_file

    query = get_dwi_file(filetype)

    assert query.pattern == expected_pattern
    assert query.description == expected_description
    assert query.needed_pipeline == expected_pipelines


@pytest.mark.parametrize(
    "filetype,expected_pattern,expected_description,expected_pipelines",
    [
        (
            "nii",
            "dwi/preprocessing/sub-*_ses-*_space-*_desc-preproc_dwi.nii*",
            "preprocessed nii files",
            "dwi-preprocessing-using-t1 or dwi-preprocessing-using-fieldmap",
        ),
        (
            "json",
            "dwi/preprocessing/sub-*_ses-*_space-*_desc-preproc_dwi.json*",
            "preprocessed json files",
            "dwi-preprocessing-using-t1 or dwi-preprocessing-using-fieldmap",
        ),
        (
            "bvec",
            "dwi/preprocessing/sub-*_ses-*_space-*_desc-preproc_dwi.bvec*",
            "preprocessed bvec files",
            "dwi-preprocessing-using-t1 or dwi-preprocessing-using-fieldmap",
        ),
        (
            "bval",
            "dwi/preprocessing/sub-*_ses-*_space-*_desc-preproc_dwi.bval*",
            "preprocessed bval files",
            "dwi-preprocessing-using-t1 or dwi-preprocessing-using-fieldmap",
        ),
    ],
)
def test_get_dwi_preprocessed_file(
    filetype: str,
    expected_pattern: str,
    expected_description: str,
    expected_pipelines: str,
):
    from clinica.utils.input_files import get_dwi_preprocessed_file

    query = get_dwi_preprocessed_file(filetype)

    assert query.pattern == expected_pattern
    assert query.description == expected_description
    assert query.needed_pipeline == expected_pipelines


def test_bids_pet_nii_empty():
    from clinica.utils.input_files import get_pet_nifti

    query = get_pet_nifti()

    assert query.pattern == str(Path("pet") / "*_pet.nii*")
    assert query.description == "PET data"


@pytest.fixture
def expected_bids_pet_query(
    tracer: Tracer, reconstruction: ReconstructionMethod
) -> QueryPattern:
    return QueryPattern(
        str(Path("pet") / f"*_trc-{tracer.value}_rec-{reconstruction.value}_pet.nii*"),
        f"PET data with {tracer.value} tracer and reconstruction method {reconstruction.value}",
        "",
    )


@pytest.mark.parametrize("tracer", Tracer)
@pytest.mark.parametrize("reconstruction", ReconstructionMethod)
def test_bids_pet_nii(
    tracer: Tracer,
    reconstruction: ReconstructionMethod,
    expected_bids_pet_query: QueryPattern,
):
    from clinica.utils.input_files import get_pet_nifti

    assert get_pet_nifti(tracer, reconstruction) == expected_bids_pet_query


@pytest.mark.parametrize("dti_measure", DTIBasedMeasure)
@pytest.mark.parametrize("space", [None, "*", "T1w"])
def test_dwi_dti_query(dti_measure, space):
    from clinica.utils.input_files import get_dwi_dti

    space = space or "*"
    query = get_dwi_dti(dti_measure, space=space)

    assert (
        query.pattern
        == f"dwi/dti_based_processing/*/*_space-{space}_{dti_measure.value}.nii.gz"
    )
    assert query.description == f"DTI-based {dti_measure.value} in space {space}."
    assert query.needed_pipeline == "dwi_dti"


def test_dwi_dti_query_error():
    from clinica.utils.input_files import get_dwi_dti

    with pytest.raises(
        ValueError,
        match="'foo' is not a valid DTIBasedMeasure",
    ):
        get_dwi_dti("foo")
