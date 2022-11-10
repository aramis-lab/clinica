from pathlib import Path

import pytest


def test_read_imaging_data(tmp_path):
    import shutil

    from clinica.iotools.converters.ukb_to_bids.ukb_utils import read_imaging_data

    path_to_zip = tmp_path / Path("alt")
    shutil.make_archive(path_to_zip, "zip", tmp_path)
    with pytest.raises(
        ValueError,
        match=f"No imaging data were found in the provided folder: {path_to_zip}, or they are not handled by Clinica. Please check your data.",
    ):
        read_imaging_data(path_to_zip)
