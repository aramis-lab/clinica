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
        detect_cross_sectional_and_longitudinal_subjects,
    )

    bids_dir = tmp_path / Path("bids")
    dirs = [
        bids_dir / Path(f"{subs[i]}") / Path(f"{directory_content[i]}")
        for i in range(0, len(subs))
    ]
    for dir in dirs:
        if not dir.is_dir():
            os.makedirs(dir)
    assert detect_cross_sectional_and_longitudinal_subjects(subs, bids_dir) == expected
