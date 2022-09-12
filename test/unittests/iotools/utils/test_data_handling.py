from pathlib import Path
from string import Template
from typing import Union

import pytest


def test_vox_to_world_space():
    """Test function `vox_to_world_space_method_1`."""

    import numpy as np
    from nibabel.nifti1 import Nifti1Image

    from clinica.iotools.utils.data_handling import vox_to_world_space_method_1

    vol = 100 + 10 * np.random.randn(5, 5, 2, 100)
    img = Nifti1Image(vol, np.eye(4))
    nii_header = img.get_header()
    coords = np.array([0.0, 0.0, 0.0])
    vox_to_world_space_method_1(coordinates_vol=coords, header=nii_header)
    assert np.array_equal(
        vox_to_world_space_method_1(coordinates_vol=coords, header=nii_header), coords
    )
