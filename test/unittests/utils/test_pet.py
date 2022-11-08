import os
from pathlib import Path

import pandas as pd
import pytest


@pytest.fixture
def psf_df() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "participant_id": ["sub-CLNC01"] * 3 + ["sub-CLNC02", "sub-CLNC03"],
            "session_id": ["ses-M000", "ses-M018"] + ["ses-M000"] * 3,
            "acq_label": ["FDG", "FDG", "AV45", "FDG", "FDG"],
            "psf_x": [8, 8, 7, 8, 8],
            "psf_y": [9, 9, 6, 9, 9],
            "psf_z": [10, 10, 5, 10, 10],
        }
    )


def test_read_psf_information_errors(tmp_path: os.PathLike, psf_df: pd.DataFrame):
    from clinica.utils.pet import read_psf_information

    with pytest.raises(
        FileNotFoundError,
        match="No such file or directory: 'foo.tsv'",
    ):
        read_psf_information(
            Path("foo.tsv"),
            ["sub-CLNC01", "sub-CLNC01"],
            ["ses-M000", "ses-M018"],
            "FDG",
        )
    psf_df.to_csv(tmp_path / "psf.tsv", sep="\t", index=False)
    with pytest.raises(
        RuntimeError,
        match=(
            "Subject sub-CLNC06 with session ses-M018 and tracer FDG "
            "that you want to proceed was not found in the TSV file containing "
            "PSF specifications"
        ),
    ):
        read_psf_information(
            tmp_path / "psf.tsv",
            ["sub-CLNC01", "sub-CLNC06"],
            ["ses-M000", "ses-M018"],
            "FDG",
        )
    psf_df_2 = pd.DataFrame(
        {
            "participant_id": ["sub-CLNC01"],
            "session_id": ["ses-M000"],
            "acq_label": ["FDG"],
            "psf_x": [10],
            "psf_y": [11],
            "psf_z": [12],
        }
    )
    duplicate_psf_df = pd.concat([psf_df, psf_df_2])
    duplicate_psf_df.to_csv(tmp_path / "duplicate_psf.tsv", sep="\t", index=False)
    with pytest.raises(
        RuntimeError,
        match=(
            "Subject sub-CLNC01 with session ses-M000 and tracer FDG "
            "that you want to proceed was found multiple times "
            "in the TSV file containing PSF specifications"
        ),
    ):
        read_psf_information(
            tmp_path / "duplicate_psf.tsv",
            ["sub-CLNC01", "sub-CLNC01"],
            ["ses-M000", "ses-M018"],
            "FDG",
        )
    psf_df["foo"] = ["bar"] * 5
    psf_df.to_csv(tmp_path / "wrong_psf.tsv", sep="\t", index=False)
    with pytest.raises(
        IOError,
        match="must contain",
    ):
        read_psf_information(
            tmp_path / "wrong_psf.tsv",
            ["sub-CLNC01", "sub-CLNC01"],
            ["ses-M000", "ses-M018"],
            "FDG",
        )
    psf_df.drop(["foo", "session_id"], axis=1).to_csv(
        tmp_path / "wrong_psf_2.tsv", sep="\t", index=False
    )
    with pytest.raises(
        IOError,
        match="must contain",
    ):
        read_psf_information(
            tmp_path / "wrong_psf_2.tsv",
            ["sub-CLNC01", "sub-CLNC01"],
            ["ses-M000", "ses-M018"],
            "FDG",
        )


def test_read_psf_information(tmp_path: os.PathLike, psf_df: pd.DataFrame):
    from clinica.utils.pet import read_psf_information

    psf_df.to_csv(tmp_path / "psf.tsv", sep="\t", index=False)
    assert read_psf_information(
        tmp_path / "psf.tsv",
        ["sub-CLNC01", "sub-CLNC01"],
        ["ses-M000", "ses-M018"],
        "FDG",
    ) == [[8, 9, 10], [8, 9, 10]]
    # Shuffle rows in dataframe and make sure results do not depend on row order
    psf_df = psf_df.sample(frac=1).reset_index(drop=True)
    psf_df.to_csv(tmp_path / "psf.tsv", sep="\t", index=False)
    assert read_psf_information(
        tmp_path / "psf.tsv",
        ["sub-CLNC01", "sub-CLNC01"],
        ["ses-M000", "ses-M018"],
        "FDG",
    ) == [[8, 9, 10], [8, 9, 10]]


@pytest.mark.parametrize(
    "label", ["pons", "cerebellumPons", "pons2", "cerebellumPons2"]
)
def test_get_suvr_mask(label: str):
    from clinica.utils.pet import get_suvr_mask

    assert Path(get_suvr_mask(label)).exists()


@pytest.mark.parametrize("label", ["foo", "bar", "pons3", "cerebelumPons2"])
def test_get_suvr_mask_error(label: str):
    from clinica.utils.pet import get_suvr_mask

    with pytest.raises(
        ValueError,
        match=f"SUVR reference region label {label} is not supported.",
    ):
        get_suvr_mask(label)
