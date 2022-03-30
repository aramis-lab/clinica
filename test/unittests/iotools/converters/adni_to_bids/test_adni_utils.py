def test_bids_id_to_loni():
    """Test function `bids_id_to_loni`."""
    from clinica.iotools.converters.adni_to_bids.adni_utils import bids_id_to_loni

    assert bids_id_to_loni("sub-ADNI000S0000") == "000_S_0000"
    assert bids_id_to_loni("sub-ADNI123S4567") == "123_S_4567"
    assert bids_id_to_loni("sub-ADNI12S4567") == "12_S_4567"
    assert bids_id_to_loni("sub-ADNI123X4567") == "123_S_4567"
    assert bids_id_to_loni("sub-ADNI123XYZ4567") == "123_S_4567"
    assert bids_id_to_loni("sub-ADNI123XYZ_TT4567") == "123_S_4567"
    assert bids_id_to_loni("sub-ADNI123XYZ12TT4567") is None
    assert bids_id_to_loni("") is None
    assert bids_id_to_loni("foo") is None
    assert bids_id_to_loni("12") is None
    assert bids_id_to_loni("123_S_4567") == "123_S_4567"
    assert bids_id_to_loni("1_XY_22") == "1_S_22"
