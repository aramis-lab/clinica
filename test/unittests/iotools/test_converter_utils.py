def test_sort_session_list():
    """Test function `sort_session_list`."""

    from clinica.iotools.converter_utils import sort_session_list

    assert sort_session_list(
        ["ses-M000", "ses-M006", "ses-M012", "ses-M024", "ses-M048", "ses-M003"]
    ) == ["ses-M000", "ses-M003", "ses-M006", "ses-M012", "ses-M024", "ses-M048"]
    assert sort_session_list(
        ["ses-M00", "ses-M06", "ses-M12", "ses-M24", "ses-M48", "ses-M03"]
    ) == ["ses-M00", "ses-M03", "ses-M06", "ses-M12", "ses-M24", "ses-M48"]
    assert sort_session_list(
        ["ses-M0", "ses-M6", "ses-M12", "ses-M24", "ses-M48", "ses-M3"]
    ) == ["ses-M0", "ses-M3", "ses-M6", "ses-M12", "ses-M24", "ses-M48"]
