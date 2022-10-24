import os
import tempfile
from pathlib import Path
from tempfile import tempdir

import pandas as pd
import pytest

from clinica.pipelines.dwi_preprocessing_using_t1.dwi_preprocessing_using_t1_utils import (
    delete_temp_dirs,
    extract_sub_ses_folder_name,
)


def test_extract_sub_ses_folder_name():
    assert (
        extract_sub_ses_folder_name(
            "/localdrive10TB/wd/dwi-preprocessing-using-t1/epi_pipeline/4336d63c8556bb56d4e9d1abc617fb3eaa3c38ea/MergeDWIs/Jacobian_image_maths_thresh_merged.nii.gz"
        )
        == "4336d63c8556bb56d4e9d1abc617fb3eaa3c38ea"
    )


@pytest.mark.parametrize(
    "checkpoint, dir_to_del, expected",
    [
        (
            Path(
                "ed42125b6abc244af649ff5eff1df41b/node_3_name/Jacobian_image_maths_thresh_merged.nii.gz"
            ),
            ["node_1_name"],
            [True, True],
        ),
        (
            Path(
                "4336d63c8556bb56d4e9d1abc617fb3eaa3c38ea/node_3_name/Jacobian_image_maths_thresh_merged.nii.gz"
            ),
            ["node_1_name"],
            [False, True],
        ),
        (
            Path(
                "4336d63c8556bb56d4e9d1abc617fb3eaa3c38ea/node_3_name/Jacobian_image_maths_thresh_merged.nii.gz"
            ),
            ["node_1_name", "node_2_name"],
            [False, False],
        ),
    ],
)
def test_delete_temp_dirs(tmp_path, checkpoint, dir_to_del, expected):
    checkpoint = tmp_path / checkpoint
    base_dir = tmp_path / Path("wd")
    dirs = [
        tmp_path
        / Path(f"wd/pipeline_name/4336d63c8556bb56d4e9d1abc617fb3eaa3c38ea/{name}")
        for name in ["node_1_name", "node_2_name"]
    ]
    for dir in dirs:
        if not dir.is_dir():
            os.makedirs(dir)
    delete_temp_dirs(checkpoint, dir_to_del, base_dir)
    assert [dir.is_dir() for dir in dirs] == expected
