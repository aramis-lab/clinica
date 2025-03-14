from pathlib import Path

import pytest

from clinica.pipelines.dwi.dti.utils import DTIBasedMeasure
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


@pytest.mark.parametrize(
    "tracer, reconstruction, expected",
    [
        (None, None, {"pattern": Path("pet/*_pet.nii*"), "description": "PET data"}),
        (
            Tracer.FDG,
            None,
            {
                "pattern": Path("pet/*_trc-18FFDG*_pet.nii*"),
                "description": "PET data with 18FFDG tracer",
            },
        ),
        (
            None,
            ReconstructionMethod.STATIC,
            {
                "pattern": Path("pet/*_rec-nacstat_pet.nii*"),
                "description": "PET data with reconstruction method nacstat",
            },
        ),
        (
            "11CPIB",
            "coregavg",
            {
                "pattern": Path("pet/*_trc-11CPIB_rec-coregavg_pet.nii*"),
                "description": "PET data with 11CPIB tracer with reconstruction method coregavg",
            },
        ),
    ],
)
def test_bids_pet_nii(tracer, reconstruction, expected):
    from clinica.utils.input_files import bids_pet_nii

    assert bids_pet_nii(tracer, reconstruction) == expected


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
