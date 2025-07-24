from cmath import nan
from pathlib import Path

import pandas as pd
import pytest


def test_read_imaging_data(tmp_path):
    import shutil

    from clinica.converters.ukb_to_bids._utils import read_imaging_data

    path_to_zip = tmp_path / "alt"
    shutil.make_archive(path_to_zip, "zip", tmp_path)
    with pytest.raises(
        ValueError,
        match=f"No imaging data were found in the provided folder: {path_to_zip}, "
        "or they are not handled by Clinica. Please check your data.",
    ):
        read_imaging_data(path_to_zip)


@pytest.mark.parametrize(
    "sidecars, expected", [([], {".nii.gz"}), (["truc.json"], {".nii.gz", ".json"})]
)
def test_get_extensions_from_sidecars_success(sidecars, expected):
    from clinica.converters.ukb_to_bids._utils import _get_extensions_from_sidecars

    assert expected == set(_get_extensions_from_sidecars(sidecars))


@pytest.mark.parametrize("sidecars", [["foo"], [".bar"]])
def test_get_extensions_from_sidecars_error(sidecars):
    from clinica.converters.ukb_to_bids._utils import _get_extensions_from_sidecars

    assert _get_extensions_from_sidecars(sidecars) == [".nii.gz"]


@pytest.mark.parametrize(
    "subject_id, source_session, age_2, age_3, expected",
    [
        ("sub1", "2", nan, nan, None),
        ("sub2", "2", 23, nan, 23),
        ("sub3", "3", nan, nan, None),
        ("sub4", "3", nan, 35, 35),
        ("sub5", "4", 0, 0, None),
    ],
)
def test_select_sessions(subject_id, source_session, age_2, age_3, expected):
    from clinica.converters.ukb_to_bids._utils import _select_sessions

    clinical_data = pd.Series(
        {
            "eid": subject_id,
            "source_sessions_number": source_session,
            "age_when_attended_assessment_centre_f21003_2_0": age_2,
            "age_when_attended_assessment_centre_f21003_3_0": age_3,
        }
    )

    assert expected == _select_sessions(clinical_data)
