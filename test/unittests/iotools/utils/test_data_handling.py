def test_vox_to_world_space():
    """Test function `vox_to_world_space_method_1`.

    .. note::
        We want to verify, that the function does run with it's normal parameters
        and that it returns what we expect it to return."""

    import numpy as np
    from nibabel.nifti1 import Nifti1Image

    from clinica.iotools.utils.data_handling import vox_to_world_space_method_1

    vol = np.ndarray(shape=(5, 5, 2, 100))
    img = Nifti1Image(vol, np.eye(4))
    coords = np.array([0.0, 0.0, 0.0])

    assert np.array_equal(
        vox_to_world_space_method_1(coordinates_vol=coords, header=img.get_header()),
        coords,
    )
