import pandas as pd
import pytest

from clinica.iotools.bids_utils import identify_modality
from clinica.iotools.converters.genfi_to_bids.genfi_to_bids_utils import compute_runs


@pytest.mark.parametrize(
    "input,expected",
    [
        ("DWI", "dwi"),
        ("T1", "T1"),
        ("T2", "T2w"),
        ("fieldmap", "fieldmap"),
        ("fmri", "rsfmri"),
    ],
)
def test_identify_modality(input, expected):
    assert identify_modality(input) == expected


def test_identify_modality_is_nan():
    import math

    assert math.isnan(identify_modality("blzflbzv"))


# test compute runs
from pandas.testing import assert_frame_equal


@pytest.mark.parametrize(
    "input, expected",
    [
        (
            pd.DataFrame(
                {
                    "source_id": ["GRN001", "GRN001", "GRN001", "C9ORF001", "C9ORF001"],
                    "source_ses_id": [1, 1, 1, 2, 2],
                    "suffix": ["a", "a", "c", "a", "b"],
                    "number_of_parts": [2, 2, 1, 1, 1],
                    "dir_num": [10, 20, 30, 10, 20],
                }
            ),
            pd.DataFrame(
                {
                    "source_id": ["C9ORF001", "C9ORF001", "GRN001", "GRN001", "GRN001"],
                    "source_ses_id": [2, 2, 1, 1, 1],
                    "suffix": ["a", "b", "a", "a", "c"],
                    "number_of_parts": [1, 1, 2, 2, 1],
                    "dir_num": [10, 20, 10, 20, 30],
                    "run_01_dir_num": [10, 20, 10, 10, 30],
                    "run": [False, False, False, True, False],
                    "run_number": [1, 1, 1, 2, 1],
                    "run_num": ["run-01", "run-01", "run-01", "run-02", "run-01"],
                }
            ),
        )
    ],
)
def test_compute_runs(input, expected):
    assert_frame_equal(compute_runs(input), expected)
