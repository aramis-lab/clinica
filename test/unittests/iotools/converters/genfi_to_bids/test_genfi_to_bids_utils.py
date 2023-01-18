import pytest

import clinica.iotools.converters.genfi_to_bids.genfi_to_bids_utils as gen


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
