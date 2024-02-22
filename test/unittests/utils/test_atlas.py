import numpy as np
import pytest
from numpy.testing import assert_array_equal

from clinica.utils.atlas import (
    AAL2,
    AICHA,
    LPBA40,
    Hammers,
    JHUDTI811mm,
    JHUTracts01mm,
    JHUTracts251mm,
    JHUTracts501mm,
    Neuromorphometrics,
)


def test_atlas_factory_error():
    from clinica.utils.atlas import atlas_factory

    with pytest.raises(
        ValueError,
        match="'foo' is not a valid AtlasName",
    ):
        atlas_factory("foo")


@pytest.mark.parametrize(
    "atlas_name,atlas",
    [
        ("AAL2", AAL2),
        ("AICHA", AICHA),
        ("Hammers", Hammers),
        ("LPBA40", LPBA40),
        ("Neuromorphometrics", Neuromorphometrics),
        ("JHUDTI81", JHUDTI811mm),
        ("JHUTract0", JHUTracts01mm),
        ("JHUTract25", JHUTracts251mm),
        ("JHUTracts50", JHUTracts501mm),
    ],
)
def test_atlas_factory(tmp_path, monkeypatch, atlas_name, atlas):
    from clinica.utils.atlas import atlas_factory

    monkeypatch.setenv("FSLDIR", str(tmp_path))
    assert isinstance(atlas_factory(atlas_name), atlas)


def test_atlas_neuromorphometrics():
    atlas = Neuromorphometrics()
    assert atlas.name == "Neuromorphometrics"
    assert (
        atlas.expected_checksum
        == "19a50136cd2f8a14357a19ad8a1dc4a2ecb6beb3fc16cb5441f4f2ebaf64a9a5"
    )
    assert atlas.atlas_filename == "atlas-Neuromorphometrics_dseg.nii.gz"
    assert atlas.roi_filename == "atlas-Neuromorphometrics_dseg.tsv"
    assert atlas.tsv_roi.exists()
    assert atlas.spatial_resolution == "1.5x1.5x1.5"
    assert atlas.atlas_folder.exists()
    assert_array_equal(atlas.get_index(), np.arange(141))
