import nibabel as nib
import numpy as np
import pytest


@pytest.mark.parametrize(
    "data_shape,expected",
    [
        ((6, 6, 6), "args"),
        ((6, 6, 290), "args -cw256"),
    ],
)
def test_check_flags(tmp_path, data_shape, expected):
    from clinica.pipelines.anatomical.freesurfer.t1.utils import _check_flags  # noqa

    nib.save(
        nib.nifti1.Nifti1Image(np.random.random(data_shape), np.eye(4)),
        tmp_path / "test.nii.gz",
    )
    assert _check_flags(tmp_path / "test.nii.gz", "args") == expected
