import os
import tempfile
from pathlib import Path
from tempfile import tempdir

import pandas as pd
import pytest

from clinica.pipelines.dwi_preprocessing_using_t1.dwi_preprocessing_using_t1_utils import (
    delete_temp_dirs,
)


@pytest.mark.parametrize(
    "checkpoint, dir_to_del, expected",
    [
        (
            Path(os.getcwd())
            / Path(
                "ed42125b6abc244af649ff5eff1df41b/node_3_name/Jacobian_image_maths_thresh_merged.nii.gz"
            ),
            ["node_1_name"],
            [True, True],
        ),
        (
            Path(os.getcwd())
            / Path(
                "4336d63c8556bb56d4e9d1abc617fb3eaa3c38ea/node_3_name/Jacobian_image_maths_thresh_merged.nii.gz"
            ),
            ["node_1_name"],
            [False, True],
        ),
        (
            Path(os.getcwd())
            / Path(
                "4336d63c8556bb56d4e9d1abc617fb3eaa3c38ea/node_3_name/Jacobian_image_maths_thresh_merged.nii.gz"
            ),
            ["node_1_name", "node_2_name"],
            [False, False],
        ),
    ],
)
def test_delete_temp_dirs(checkpoint, dir_to_del, expected):

    base_path = os.getcwd()
    dir_1 = Path(base_path) / Path(
        "wd/pipeline_name/4336d63c8556bb56d4e9d1abc617fb3eaa3c38ea/node_1_name"
    )
    dir_2 = Path(base_path) / Path(
        "wd/pipeline_name/4336d63c8556bb56d4e9d1abc617fb3eaa3c38ea/node_2_name"
    )
    if not Path(dir_1).is_dir():
        os.makedirs(dir_1)
    if not Path(dir_2).is_dir():
        os.makedirs(dir_2)
    base_dir = Path(base_path) / Path("wd")
    delete_temp_dirs(checkpoint, dir_to_del, base_dir)
    assert Path(dir_1).is_dir() == expected[0]
    assert Path(dir_2).is_dir() == expected[1]
