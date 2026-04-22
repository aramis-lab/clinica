from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from pandas.testing import assert_frame_equal, assert_series_equal

from clinica.converters.oasis3_to_bids._utils import (
    _find_csv_with_filename,
    _read_clinical_data,
)


def _get_mri_data() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "OASISID": ["OAS30001", "OAS30002"],
            "label": ["OAS30001_MR_d0000", "OAS30002_MR_d0101"],
            "Manufacturer": ["Siemens", "Siemens"],
            "ManufacturersModelName": ["Biograph_mMR", "Biograph_mMR"],
            "MagneticFieldStrength": [3, 3],
        }
    )


def _get_imaging_data() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "Subject": ["OAS30001", "OAS30001"],
            "Date": ["d0000", "d0200"],
            "path": ["OAS30001_MR_d0000", "OAS30001_MR_d0200"],
        }
    )


def _get_demo_data() -> pd.DataFrame:
    return pd.DataFrame(
        data={
            "OASISID": ["OAS30001", "OAS30002"],
            "AgeatEntry": [62, 45],
            "GENDER": [2, 1],
        }
    )


def _build_clinical_data(tmp_path: Path) -> Path:
    clinical_data_directory = tmp_path / "clinical_data"
    clinical_data_directory.mkdir()

    udsb1 = pd.DataFrame(
        data={
            "OASISID": ["OAS30001", "OAS30001", "OAS30002"],
            "days_to_visit": ["0000", "0200", "0101"],
            "age at visit": [62, 62.55, 45.67],
            "WEIGHT": [152, 152, 184.2],
            "MMSE": [28, 28, 29],
            "CDRTOT": [0, 0, 1],
            "commun": [0, 0, 0],
            "dx1": [1, 1, 1],
            "homehobb": [0, 0, 0],
            "judgment": [0, 0, 0],
            "orient": [0, 0, 0],
            "perscare": [0, 0, 0],
            "CDRSUM": [0, 0, 1],
        }
    )

    udsb1.to_csv(
        clinical_data_directory / "OASIS3_UDSb1_physical_eval.csv", index=False
    )

    udsb4 = pd.DataFrame(
        data={
            "OASISID": ["OAS30001", "OAS30001", "OAS30002"],
            "days_to_visit": ["0000", "0200", "0101"],
            "age at visit": [62, 62.55, 45.67],
            "memory": [0, 0, 1],
        }
    )
    udsb4.to_csv(clinical_data_directory / "OASIS3_UDSb4_cdr.csv", index=False)

    _get_demo_data().to_csv(
        clinical_data_directory / "OASIS3_demographics.csv", index=False
    )

    return clinical_data_directory


def test_find_csv_with_filename_success(tmp_path):
    expected = pd.DataFrame({"test_data": ["foo", "bar", "foobar"]})
    expected.to_csv(tmp_path / "foo.csv", index=False)
    assert_frame_equal(_find_csv_with_filename(tmp_path, filename="foo"), expected)


def test_find_csv_with_filename_no_file(tmp_path):
    with pytest.raises(FileNotFoundError, match=f"No CSV files found"):
        _find_csv_with_filename(tmp_path, "foo")


def test_find_csv_with_filename_more_files(tmp_path):
    (tmp_path / "foo.csv").touch()
    (tmp_path / "dupe_foo.csv").touch()
    with pytest.raises(FileNotFoundError, match="More than one file"):
        _find_csv_with_filename(tmp_path, filename="foo")


def test_read_mri_data(tmp_path):
    from clinica.converters.oasis3_to_bids._utils import _read_mri_data

    pd.DataFrame().to_csv(tmp_path / "OASIS3_MR_json.csv")
    _read_mri_data(tmp_path)


def test_find_imaging_data(tmp_path):
    from clinica.converters.oasis3_to_bids._utils import _find_imaging_data

    (tmp_path / "images").mkdir()
    (tmp_path / "images" / "test.nii.gz").touch()
    assert "images/test.nii.gz" == str(next(_find_imaging_data(tmp_path)))


@pytest.mark.parametrize(
    "input, expected",
    [
        ("sub-01_ses-02_t1.nii.gz", "run-01"),
        ("imaging_data/sub-01_ses-02_run-12_t1.nii.gz", "run-12"),
    ],
)
def test_identify_runs(input, expected):
    from clinica.converters.oasis3_to_bids._utils import _identify_runs

    assert expected == _identify_runs(Path(input))


def test_read_clinical_data(tmp_path):
    result = _read_clinical_data(_build_clinical_data(tmp_path))
    expected = pd.DataFrame(
        data={
            "Subject": ["OAS30001", "OAS30001", "OAS30002"],
            "days_to_visit": [0, 200, 101],
            "age at visit": [62, 62.55, 45.67],
            "WEIGHT": [152, 152, 184.2],
            "memory": [0, 0, 1],
            "Date": ["d0000", "d0200", "d0101"],
            "mmse": [28, 28, 29],
            "cdr": [0, 0, 1],
            "commun": [0, 0, 0],
            "dx1": [1, 1, 1],
            "homehobb": [0, 0, 0],
            "judgment": [0, 0, 0],
            "orient": [0, 0, 0],
            "perscare": [0, 0, 0],
            "sumbox": [0, 0, 1],
        }
    )
    assert_frame_equal(expected, result, check_like=True)


def test_read_demo_data(tmp_path):
    from clinica.converters.oasis3_to_bids._utils import _read_demo_data

    expected = pd.DataFrame(
        data={
            "Subject": ["OAS30001", "OAS30002"],
            "AgeatEntry": [62, 45],
            "sex": ["F", "M"],
        }
    )
    result = _read_demo_data(_build_clinical_data(tmp_path))
    assert_frame_equal(expected, result, check_exact=False)


def test_map_modality_to_bids():
    from clinica.converters.oasis3_to_bids._utils import _map_supported_modality_to_bids

    input = pd.DataFrame(
        {
            "modality": [
                "dwi_MR",
                "truc",
                "T1w_MR",
                "T2star_MR",
                "FLAIR_MR",
                "pet_FDG",
                "pet_PIB",
                "pet_AV45",
                "bold_MR",
                "pet_AV1451",
            ]
        }
    )
    expected = pd.DataFrame(
        {
            "modality": [
                "dwi_MR",
                "truc",
                "T1w_MR",
                "T2star_MR",
                "FLAIR_MR",
                "pet_FDG",
                "pet_PIB",
                "pet_AV45",
                "bold_MR",
                "pet_AV1451",
            ],
            "datatype": [
                "dwi",
                np.nan,
                "anat",
                "anat",
                "anat",
                "pet",
                "pet",
                "pet",
                "func",
                "pet",
            ],
            "suffix": [
                "dwi",
                np.nan,
                "T1w",
                "T2starw",
                "FLAIR",
                "pet",
                "pet",
                "pet",
                "bold",
                "pet",
            ],
            "trc_label": [
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                "18FFDG",
                "11CPIB",
                "18FAV45",
                np.nan,
                "18FAV1451",
            ],
        }
    )

    assert_frame_equal(
        _map_supported_modality_to_bids(input), expected, check_like=True
    )


def test_filter_to_imaging_subjects():
    from clinica.converters.oasis3_to_bids._utils import (
        _filter_to_imaging_subjects,
    )

    demo_data = _get_demo_data().rename({"OASISID": "Subject"}, axis=1)
    result = _filter_to_imaging_subjects(demo_data, _get_imaging_data()["Subject"])
    expected = pd.DataFrame(
        data={
            "Subject": ["OAS30001"],
            "AgeatEntry": [62],
            "GENDER": [2],
        }
    )
    assert_frame_equal(result, expected, check_like=True)


def test_get_baseline_visit_clinical(tmp_path):
    from clinica.converters.oasis3_to_bids._utils import (
        _get_baseline_visit_clinical,
    )

    clinical_data = _read_clinical_data(_build_clinical_data(tmp_path)).rename(
        {"OASISID": "Subject"}, axis=1
    )
    expected = pd.DataFrame(
        {
            "Subject": ["OAS30001"],
            "days_to_visit": [0],
            "age at visit": [62.0],
            "WEIGHT": [152.0],
            "memory": [0],
            "Date": ["d0000"],
            "mmse": [28],
            "cdr": [0],
            "commun": [0],
            "dx1": [1],
            "homehobb": [0],
            "judgment": [0],
            "orient": [0],
            "perscare": [0],
            "sumbox": [0],
        }
    )
    result = _get_baseline_visit_clinical(clinical_data, _get_imaging_data()["Subject"])
    assert_frame_equal(result, expected, check_like=True)


def test_add_age_at_entry_and_age_at_scan():
    from clinica.converters.oasis3_to_bids._utils import (
        _add_age_at_entry_and_age_at_scan,
    )

    demo_data = _get_demo_data().rename({"OASISID": "Subject"}, axis=1)
    result = _add_age_at_entry_and_age_at_scan(_get_imaging_data(), demo_data)
    expected = pd.DataFrame(
        {
            "Subject": ["OAS30001", "OAS30001"],
            "Date": ["d0000", "d0200"],
            "AgeatEntry": [62, 62],
            "age": [62, 62.55],
            "path": ["OAS30001_MR_d0000", "OAS30001_MR_d0200"],
        }
    )
    assert_frame_equal(result, expected, check_like=True)


def test_add_mri_scanner_metadata():
    from clinica.converters.oasis3_to_bids._utils import _add_mri_scanner_metadata

    result = _add_mri_scanner_metadata(_get_imaging_data(), _get_mri_data())
    expected = pd.DataFrame(
        {
            "Subject": ["OAS30001", "OAS30001"],
            "Date": ["d0000", "d0200"],
            "path": ["OAS30001_MR_d0000", "OAS30001_MR_d0200"],
            "label": ["OAS30001_MR_d0000", np.nan],
            "Manufacturer": ["Siemens", np.nan],
            "ManufacturersModelName": ["Biograph_mMR", np.nan],
            "MagneticFieldStrength": [3, np.nan],
        }
    )
    assert_frame_equal(result, expected, check_like=True)


def test_merge_baseline_data():
    from clinica.converters.oasis3_to_bids._utils import _merge_baseline_data

    input1 = pd.DataFrame(
        data={
            "Subject": ["OAS30001", "OAS30002"],
            "AgeatEntry": [62, 45],
        }
    )
    input2 = pd.DataFrame(
        data={
            "Subject": ["OAS30001", "OAS30004"],
            "sex": ["F", "M"],
        }
    )
    result = _merge_baseline_data(input1, input2)
    expected = pd.DataFrame(
        data={
            "Subject": ["OAS30001"],
            "AgeatEntry": [62],
            "sex": ["F"],
            "participant_id": ["sub-OAS30001"],
        }
    )
    assert_frame_equal(result, expected, check_like=True)


def test_days_to_session_bin():
    from clinica.converters.oasis3_to_bids._utils import _days_to_session_bin

    assert_series_equal(
        _days_to_session_bin(pd.Series(data=[0, 200, 300, 360])),
        pd.Series([0, 6, 12, 12]),
    )


def test_compute_session_bins():
    from clinica.converters.oasis3_to_bids._utils import _compute_session_bins

    result = _compute_session_bins(_get_imaging_data())
    expected = pd.DataFrame(
        {
            "Subject": ["OAS30001", "OAS30001"],
            "Date": ["d0000", "d0200"],
            "path": ["OAS30001_MR_d0000", "OAS30001_MR_d0200"],
            "session": [0, 6],
            "ses": ["ses-M000", "ses-M006"],
        }
    )
    assert_frame_equal(result, expected, check_like=True)


def test_merge_clinical_scores(tmp_path):
    from clinica.converters.oasis3_to_bids._utils import (
        _compute_session_bins,
        _merge_clinical_scores,
    )

    result = _merge_clinical_scores(
        _compute_session_bins(_get_imaging_data()),
        _read_clinical_data(_build_clinical_data(tmp_path)),
    )
    expected = pd.DataFrame(
        data={
            "Subject": ["OAS30001", "OAS30001"],
            "Date": ["d0000", "d0200"],
            "path": ["OAS30001_MR_d0000", "OAS30001_MR_d0200"],
            "session": [0, 6],
            "ses": ["ses-M000", "ses-M006"],
            "memory": [0, 0],
            "mmse": [28, 28],
            "cdr": [0, 0],
            "commun": [0, 0],
            "dx1": [1, 1],
            "homehobb": [0, 0],
            "judgment": [0, 0],
            "orient": [0, 0],
            "perscare": [0, 0],
            "sumbox": [0, 0],
        }
    )
    assert_frame_equal(result, expected, check_like=True)


def test_filter_imaging_data_by_subjects(tmp_path):
    from clinica.converters.oasis3_to_bids._utils import filter_imaging_data_by_subjects

    imaging_data = pd.concat(
        [
            _get_imaging_data(),
            pd.DataFrame(
                {
                    "Subject": ["OAS30002"],
                    "Date": ["d0000"],
                    "path": ["OAS30002_MR_d0000"],
                }
            ),
        ]
    )
    with open(tmp_path / "subjects.txt", "w") as f:
        f.write("OAS30001\nfoo\n")

    output = filter_imaging_data_by_subjects(imaging_data, tmp_path / "subjects.txt")
    assert len(output) == 2
    assert output.Subject.unique()[0] == "OAS30001"
