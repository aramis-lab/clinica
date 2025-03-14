import math

import numpy as np
import pandas as pd
import pytest
from pandas.testing import assert_frame_equal


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
def test_identify_modality(input: str, expected: str):
    from clinica.converters._utils import identify_modality

    assert identify_modality(input) == expected


def test_identify_modality_is_nan():
    from clinica.converters._utils import identify_modality

    assert math.isnan(identify_modality("blzflbzv"))


@pytest.fixture
def input_df_compute_runs() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "source_id": ["GRN001", "GRN001", "GRN001", "C9ORF001", "C9ORF001"],
            "source_ses_id": [1, 1, 1, 2, 2],
            "suffix": ["a", "a", "c", "a", "b"],
            "number_of_parts": [2, 2, 1, 1, 1],
            "dir_num": [10, 20, 30, 10, 20],
        }
    )


def test_compute_runs(input_df_compute_runs: pd.DataFrame):
    from clinica.converters.genfi_to_bids._utils import _compute_run_numbers_from_parts

    expected = pd.DataFrame(
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
    )
    assert_frame_equal(_compute_run_numbers_from_parts(input_df_compute_runs), expected)


@pytest.fixture
def input_df_compute_philips_parts() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "source_id": ["sub-01"] * 5 + ["sub-02"] * 3,
            "source_ses_id": [1, 1, 1, 2, 2, 1, 2, 2],
            "suffix": ["dwi"] * 8,
            "dir_num": [10, 20, 30, 40, 50, 10, 20, 30],
        }
    )


def test_compute_philips_parts(input_df_compute_philips_parts: pd.DataFrame):
    from clinica.converters.genfi_to_bids._utils import _compute_philips_parts

    expected = pd.DataFrame(
        {
            "source_id": ["sub-01"] * 5 + ["sub-02"] * 3,
            "source_ses_id": [1, 1, 1, 2, 2, 1, 2, 2],
            "suffix": ["dwi"] * 8,
            "dir_num": [10, 20, 30, 40, 50, 10, 20, 30],
            "part_01_dir_num": [10, 10, 10, 40, 40, 10, 20, 20],
            "run": [False, True, True, False, True, False, False, True],
            "part_number": [1, 2, 3, 1, 2, 1, 1, 2],
            "number_of_parts": [3, 3, 3, 2, 2, 1, 2, 2],
        }
    )
    assert_frame_equal(_compute_philips_parts(input_df_compute_philips_parts), expected)


@pytest.mark.parametrize(
    "input, expected",
    [
        (
            [False, True, False, True, False, False, False, True, False, True, True],
            [1, 2, 1, 2, 1, 1, 1, 2, 1, 2, 3],
        ),
        ([True, True, True], [1, 2, 3]),
        ([False, False, False], [1, 1, 1]),
    ],
)
def test_compute_scan_sequence_numbers(input: list[bool], expected: list[int]):
    from clinica.converters.genfi_to_bids._utils import _compute_scan_sequence_numbers

    assert _compute_scan_sequence_numbers(input) == expected


def test_compute_scan_sequence_numbers_error():
    from clinica.converters.genfi_to_bids._utils import _compute_scan_sequence_numbers

    with pytest.raises(ValueError):
        _compute_scan_sequence_numbers([])


def test_drop_duplicate_line_with_nans_error():
    from clinica.converters.genfi_to_bids._utils import _drop_duplicate_line_with_nans

    with pytest.raises(
        ValueError,
        match="Column participant_id was not found in the participants tsv while it is required by BIDS specifications.",
    ):
        _drop_duplicate_line_with_nans(pd.DataFrame())


def test_drop_duplicate_line_with_nans():
    from clinica.converters.genfi_to_bids._utils import _drop_duplicate_line_with_nans

    assert_frame_equal(
        _drop_duplicate_line_with_nans(
            pd.DataFrame(
                {
                    "participant_id": ["sub-1", "sub-1", "sub-2", "sub-3"],
                    "metadata": ["1", np.nan, np.nan, "3"],
                }
            )
        ).reset_index(drop=True),
        pd.DataFrame(
            {
                "participant_id": ["sub-1", "sub-2", "sub-3"],
                "metadata": ["1", np.nan, "3"],
            }
        ),
    )
