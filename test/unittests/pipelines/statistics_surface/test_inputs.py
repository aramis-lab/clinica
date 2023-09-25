import os
from pathlib import Path

import pandas as pd
import pytest

CURRENT_DIR = os.path.dirname(os.path.realpath(__file__))


def test_read_and_check_tsv_file_filenotfound_error(tmpdir):
    from clinica.pipelines.statistics_surface.surfstat._utils import (
        read_and_check_tsv_file,
    )

    with pytest.raises(FileNotFoundError, match="File foo.tsv does not exist"):
        read_and_check_tsv_file(Path("foo.tsv"))


@pytest.mark.parametrize(
    "columns", [["foo"], ["foo", "bar"], ["participant_id", "bar"]]
)
def test_read_and_check_tsv_file_data_errors(tmpdir, columns):
    from clinica.pipelines.statistics_surface.surfstat._utils import (
        read_and_check_tsv_file,
    )

    df = pd.DataFrame(columns=columns)
    df.to_csv(tmpdir / "foo.tsv", sep="\t", index=False)
    with pytest.raises(
        ValueError,
        match=r"The TSV data should have at least two columns: participant_id and session_id",
    ):
        read_and_check_tsv_file(tmpdir / "foo.tsv")


def test_read_and_check_tsv_file():
    from clinica.pipelines.statistics_surface.surfstat._utils import (
        read_and_check_tsv_file,
    )

    df = read_and_check_tsv_file(Path(CURRENT_DIR) / "data/subjects.tsv")
    assert len(df) == 7
    assert set(df.columns) == {"group", "age", "sex"}
