import pandas as pd
import pytest
from pandas.testing import assert_frame_equal

from clinica.converters.oasis3_to_bids._utils import _find_csv_with_filename


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
