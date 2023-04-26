import pandas as pd
import pytest

from clinica.iotools.bids_utils import identify_modality


@pytest.mark.parametrize(
    "input,expected",
    [
        ("DWI", "dwi"),
        ("T1", "T1"),
        ("T2", "T2w"),
        ("fieldmap", "fieldmap"),
        ("fmri", "rsfmri"),
        # ("blzflbzv", pd.NA),
    ],
)
def test_identify_modality(input, expected):
    assert identify_modality(input) == expected


def test_identify_modality_is_nan():
    import math

    assert math.isnan(identify_modality("blzflbzv"))
