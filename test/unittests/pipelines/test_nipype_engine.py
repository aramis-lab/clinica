import os
from pathlib import Path

import pytest


@pytest.mark.parametrize(
    "subs, directory_content, expected",
    [
        (["sub-01", "sub-02"], ["ses-M001", "ses-M001"], ([], ["sub-01", "sub-02"])),
        (["sub-01", "sub-02"], ["anat", "ses-M001"], (["sub-01"], ["sub-02"])),
        (["sub-01", "sub-02"], ["anat", "anat"], (["sub-01", "sub-02"], [])),
    ],
)
def test_detect_cross_sectional_and_longitudinal_subjects(
    tmp_path, subs, directory_content, expected
):
    from clinica.pipelines.engine import (
        _detect_cross_sectional_and_longitudinal_subjects,
    )

    bids_dir = tmp_path / "bids"
    directories = [
        bids_dir / f"{subject}" / f"{directory_content[i]}"
        for i, subject in enumerate(subs)
    ]
    for directory in directories:
        if not directory.is_dir():
            directory.mkdir(parents=True)

    assert _detect_cross_sectional_and_longitudinal_subjects(subs, bids_dir) == expected


def _build_cross_sectional_bids(folder: Path):
    folder.mkdir(exist_ok=True, parents=True)
    for filename in ("dataset_description.json", ".bidsignore", "foo.txt"):
        (folder / filename).touch()
    for subject in ("sub-01", "sub-02"):
        (folder / subject / "anat").mkdir(parents=True, exist_ok=True)
        for extension in (".nii.gz", ".json"):
            (folder / subject / "anat" / f"{subject}_T1w{extension}").touch()


def test_convert_cross_sectional(tmp_path):
    from clinica.pipelines.engine import _convert_cross_sectional

    bids_in = tmp_path / "bids_in"
    _build_cross_sectional_bids(bids_in)
    bids_out = tmp_path / "bids_out"
    bids_out.mkdir()

    _convert_cross_sectional(bids_in, bids_out, ("sub-01", "sub-02"), ())

    assert set([f.name for f in bids_out.iterdir()]) == {
        "dataset_description.json",
        ".bidsignore",
        "sub-01",
        "sub-02",
    }
    for subject in ("sub-01", "sub-02"):
        for extension in (".nii.gz", ".json"):
            assert (
                bids_out
                / subject
                / "ses-M000"
                / "anat"
                / f"{subject}_ses-M000_T1w{extension}"
            ).exists()
