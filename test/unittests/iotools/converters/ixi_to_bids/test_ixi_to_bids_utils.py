from pathlib import Path
from typing import List

import pandas as pd
import pytest

from clinica.iotools.converters.ixi_to_bids.ixi_to_bids_utils import (
    get_subjects_list_from_data,
)


def test_get_subjects_list_from_data(tmp_path):
    for filename in ["IXI1", "IXI123", "foo"]:
        Path(f"{tmp_path}/{filename}_T1w.nii.gz").touch()
    assert get_subjects_list_from_data(tmp_path) == ["IXI123"]
