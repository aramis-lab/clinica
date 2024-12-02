from pathlib import Path

import pytest

from clinica.utils.dwi import DTIBasedMeasure
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


def test_bids_pet_nii_empty():
    from clinica.utils.input_files import bids_pet_nii

    assert bids_pet_nii() == {
        "pattern": Path("pet") / "*_pet.nii*",
        "description": "PET data",
    }


@pytest.fixture
def expected_bids_pet_query(tracer, reconstruction):
    return {
        "pattern": Path("pet")
        / f"*_trc-{tracer.value}_rec-{reconstruction.value}_pet.nii*",
        "description": f"PET data with {tracer.value} tracer and reconstruction method {reconstruction.value}",
    }


@pytest.mark.parametrize("tracer", Tracer)
@pytest.mark.parametrize("reconstruction", ReconstructionMethod)
def test_bids_pet_nii(tracer, reconstruction, expected_bids_pet_query):
    from clinica.utils.input_files import bids_pet_nii

    assert bids_pet_nii(tracer, reconstruction) == expected_bids_pet_query


@pytest.mark.parametrize("dti_measure", DTIBasedMeasure)
@pytest.mark.parametrize("space", [None, "*", "T1w"])
def test_dwi_dti_query(dti_measure, space):
    from clinica.utils.input_files import dwi_dti

    space = space or "*"
    assert dwi_dti(dti_measure, space=space) == {
        "pattern": f"dwi/dti_based_processing/*/*_space-{space}_{dti_measure.value}.nii.gz",
        "description": f"DTI-based {dti_measure.value} in space {space}.",
        "needed_pipeline": "dwi_dti",
    }


def test_dwi_dti_query_error():
    from clinica.utils.input_files import dwi_dti

    with pytest.raises(
        ValueError,
        match="'foo' is not a valid DTIBasedMeasure",
    ):
        dwi_dti("foo")


@pytest.mark.parametrize(
    "input_parameters,expected_description,expected_pattern",
    [
        (
            {},
            (
                "Cropped PET nifti image affinely registered to the MNI152NLin2009cSym "
                "template resulting from the pet-linear pipeline"
            ),
            "pet_linear/*_pet_space-MNI152NLin2009cSym_desc-Crop*_pet.nii.gz",
        ),
        (
            {
                "acq_label": "18FFDG",
            },
            (
                "Cropped PET nifti image obtained with tracer 18FFDG affinely registered to the "
                "MNI152NLin2009cSym template resulting from the pet-linear pipeline"
            ),
            "pet_linear/*_trc-18FFDG_pet_space-MNI152NLin2009cSym_desc-Crop*_pet.nii.gz",
        ),
        (
            {
                "acq_label": "18FAV45",
                "suvr_reference_region": "pons",
                "uncropped_image": True,
                "resolution": 2,
            },
            (
                "PET nifti image of resolution 2x2x2 obtained with tracer 18FAV45 for SUVR region "
                "pons affinely registered to the MNI152NLin2009cSym template resulting from the pet-linear pipeline"
            ),
            "pet_linear/*_trc-18FAV45_pet_space-MNI152NLin2009cSym*_res-2x2x2_suvr-pons_pet.nii.gz",
        ),
        (
            {
                "suvr_reference_region": "pons2",
                "resolution": 2,
                "space": "T1w",
            },
            (
                "PET nifti image of resolution 2x2x2 for SUVR region pons2 affinely "
                "registered to the associated T1w image resulting from the pet-linear pipeline"
            ),
            "pet_linear/*_pet_space-T1w*_res-2x2x2_suvr-pons2_pet.nii.gz",
        ),
        (
            {
                "acq_label": "18FFDG",
                "space": "T1w",
            },
            (
                "PET nifti image obtained with tracer 18FFDG affinely registered to the "
                "associated T1w image resulting from the pet-linear pipeline"
            ),
            "pet_linear/*_trc-18FFDG_pet_space-T1w*_pet.nii.gz",
        ),
    ],
)
def test_pet_linear_nii(input_parameters, expected_description, expected_pattern):
    from clinica.utils.input_files import pet_linear_nii

    query = pet_linear_nii(**input_parameters)

    assert query["description"] == expected_description
    assert query["needed_pipeline"] == "pet-linear"
    assert str(query["pattern"]) == expected_pattern
