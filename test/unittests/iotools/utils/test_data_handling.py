import pytest


def test_vox_to_world_space_method_1():
    """Test function `vox_to_world_space_method_1`."""

    import numpy as np
    from nibabel.nifti1 import Nifti1Header

    from clinica.iotools.utils.data_handling import vox_to_world_space_method_1

    head = Nifti1Header()
    assert not head["qform_code"] > 0
    assert head["sform_code"] == 0
    assert np.all(head["pixdim"] == 1)
    center_coordinates = np.array([0.5] * 3)
    np.testing.assert_array_equal(
        vox_to_world_space_method_1(coordinates_vol=center_coordinates, header=head),
        center_coordinates,
    )
    head["pixdim"][1] = 2
    np.testing.assert_array_equal(
        vox_to_world_space_method_1(coordinates_vol=center_coordinates, header=head),
        np.array([1.0, 0.5, 0.5]),
    )


def test_validate_output_tsv_path(tmp_path):
    from clinica.iotools.utils.data_handling import _validate_output_tsv_path

    (tmp_path / "foo.tsv").touch()

    assert _validate_output_tsv_path(tmp_path) == tmp_path / "merge.tsv"
    assert _validate_output_tsv_path(tmp_path / "foo.tsv") == tmp_path / "foo.tsv"
    assert (
        _validate_output_tsv_path(tmp_path / "bids" / "foo.tsv")
        == tmp_path / "bids" / "foo.tsv"
    )
    assert (tmp_path / "bids").exists()
    assert not (tmp_path / "bids" / "foo.tsv").exists()


def test_validate_output_tsv_path_error(tmp_path):
    from clinica.iotools.utils.data_handling import _validate_output_tsv_path

    (tmp_path / "bar.txt").touch()

    with pytest.raises(
        TypeError,
        match="Output path extension must be tsv.",
    ):
        _validate_output_tsv_path(tmp_path / "bar.txt")
