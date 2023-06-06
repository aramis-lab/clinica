def test_get_group_message_empty(tmp_path):
    from clinica.utils.ux import _get_group_message

    assert _get_group_message(tmp_path) == "No group was found in CAPS directory"

    (tmp_path / "foo.nii.gz").touch()

    assert _get_group_message(tmp_path) == "No group was found in CAPS directory"


def test_get_group_message(tmp_path):
    from clinica.utils.ux import _get_group_message

    (tmp_path / "groups").mkdir()
    for group in ("group-AD", "group-foo", "group-foobar66"):
        (tmp_path / "groups" / group).mkdir()

    assert (
        _get_group_message(tmp_path)
        == "Groups that exist in your CAPS directory are AD, foo, foobar66."
    )
