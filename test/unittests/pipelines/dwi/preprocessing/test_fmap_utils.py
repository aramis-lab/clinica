import math
from pathlib import Path

import nibabel as nib
import numpy as np
from numpy.testing import assert_array_almost_equal


def test_get_output_file(tmp_path):
    from clinica.pipelines.dwi.preprocessing.fmap.utils import _get_output_file  # noqa

    assert _get_output_file(tmp_path / "tmp" / "foo.nii.gz", "bar") == "foo_bar.nii.gz"


def write_input_image(output_dir: Path, filename: str = "foo.nii.gz"):
    img_data = np.zeros((5, 5, 5, 8))
    img_data[2:4, 2:4, 2:4, 0:2] = 1.0
    img = nib.Nifti1Image(img_data, affine=np.eye(4))
    nib.save(img, output_dir / filename)


def test_convert_phase_difference_to_hertz(tmp_path):
    from clinica.pipelines.dwi.preprocessing.fmap.utils import (
        convert_phase_difference_to_hertz,
    )

    (tmp_path / "tmp").mkdir()
    write_input_image(tmp_path)
    output = convert_phase_difference_to_hertz(
        phase_diff_filename=tmp_path / "foo.nii.gz",
        delta_echo_time=2.0,
        working_dir=tmp_path / "tmp",
    )

    assert output == tmp_path / "tmp" / "foo_radsec.nii.gz"

    output_img = nib.load(output)
    input_img = nib.load(tmp_path / "foo.nii.gz")
    assert_array_almost_equal(
        output_img.get_fdata(),
        input_img.get_fdata().astype(np.float32) * (1.0 / (4.0 * math.pi)),
    )


def test_demean_image(tmp_path):
    from clinica.pipelines.dwi.preprocessing.fmap.utils import demean_image

    (tmp_path / "tmp").mkdir()
    write_input_image(tmp_path)

    output = demean_image(
        input_image=tmp_path / "foo.nii.gz", working_dir=tmp_path / "tmp"
    )

    assert output == tmp_path / "tmp" / "foo_demean.nii.gz"

    output_img = nib.load(output)
    input_img = nib.load(tmp_path / "foo.nii.gz")
    data = input_img.get_fdata().astype(np.float32)
    data = data - np.median(data.reshape(-1))
    assert_array_almost_equal(output_img.get_fdata(), data)


def test_convert_phase_difference_to_rads(tmp_path):
    from clinica.pipelines.dwi.preprocessing.fmap.utils import (
        convert_phase_difference_to_rads,
    )

    (tmp_path / "tmp").mkdir()
    write_input_image(tmp_path)

    output = convert_phase_difference_to_rads(
        phase_diff_filename=tmp_path / "foo.nii.gz", working_dir=tmp_path / "tmp"
    )

    assert output == tmp_path / "tmp" / "foo_rads.nii.gz"
