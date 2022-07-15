import os
from pathlib import Path

import pandas as pd
import pytest

CURRENT_DIR = os.path.dirname(os.path.realpath(__file__))


def test_read_and_check_tsv_file(tmpdir):
    from clinica.pipelines.statistics_surface._inputs import _read_and_check_tsv_file

    with pytest.raises(FileNotFoundError, match="File foo.tsv does not exist"):
        _read_and_check_tsv_file(Path("foo.tsv"))
    df = pd.DataFrame(columns=["foo"])
    df.to_csv(tmpdir / "foo.tsv", sep="\t", index=False)
    with pytest.raises(
        ValueError, match=r"The TSV data in .*foo.tsv should have at least 2 columns."
    ):
        _read_and_check_tsv_file(tmpdir / "foo.tsv")
    df = pd.DataFrame(columns=["foo", "bar"])
    df.to_csv(tmpdir / "foo.tsv", sep="\t", index=False)
    with pytest.raises(
        ValueError,
        match=r"The first column in .*foo.tsv should always be participant_id.",
    ):
        _read_and_check_tsv_file(tmpdir / "foo.tsv")
    df = pd.DataFrame(columns=["participant_id", "bar"])
    df.to_csv(tmpdir / "foo.tsv", sep="\t", index=False)
    with pytest.raises(
        ValueError, match=r"The second column in .*foo.tsv should always be session_id."
    ):
        _read_and_check_tsv_file(tmpdir / "foo.tsv")
    df = _read_and_check_tsv_file(Path(CURRENT_DIR) / "data/subjects.tsv")
    assert len(df) == 7
    assert set(df.columns) == {"participant_id", "session_id", "group", "age", "sex"}
