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


def test_write_row_in_scans_tsv_file(tmp_path):
    from clinica.converters.ukb_to_bids._utils import write_row_in_scans_tsv_file

    row = pd.Series(
        {
            "participant_id": "sub-0001",
            "sessions": "ses-M000",
            "filename": "sub-0001_ses-M000_T1w.nii.gz",
            "modality": "T1w",
        }
    )

    target_dir = tmp_path / "BIDS" / "sub-0001" / "ses-M000"
    target_dir.mkdir(parents=True)

    write_row_in_scans_tsv_file(row, tmp_path / "BIDS")

    scans_tsv = target_dir / "sub-0001_ses-M000_scans.tsv"
    assert scans_tsv.exists()

    content = scans_tsv.read_text().strip().splitlines()

    columns_names = content[0].split("\t")
    columns_items = content[1].split("\t")

    assert columns_names == ["filename", "modality"]
    assert columns_items == ["sub-0001_ses-M000_T1w.nii.gz", "T1w"]


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
