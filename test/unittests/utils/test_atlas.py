import pytest

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
def test_atlas_factory(atlas_name, atlas):
    from clinica.utils.atlas import atlas_factory

    assert isinstance(atlas_factory(atlas_name)(), atlas)
