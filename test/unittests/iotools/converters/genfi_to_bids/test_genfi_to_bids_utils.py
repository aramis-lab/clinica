import pytest

from clinica.iotools.bids_utils import identify_modality


@pytest.mark.parametrize(
    "input,expected",
    [
        ("DWI", "dwi"),
        ("T1", "T1"),
        ("T2", "T2w"),
        ("fieldmap", "fieldmap"),
        ("fmri", "rsfmri"),
        ("blzflbzv", None),
    ],
)
def test_identify_modality(input, expected):
    assert gen.identify_modality(input) == expected
