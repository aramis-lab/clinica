import pytest


@pytest.mark.parametrize(
    "label",
    ["A", "a1", "123", "0xY2", "mygroup"],
)
def test_check_group_label(label):
    from clinica.utils.group import check_group_label

    assert check_group_label(label) is None


@pytest.mark.parametrize(
    "label",
    [
        "",
        "$",
        "@group1",
        "group1+group2",
        "my group",
        " ",
        "#12",
        "my_group",
        "my_group2",
    ],
)
def test_check_group_label_errors(label):
    from clinica.utils.group import check_group_label

    with pytest.raises(
        ValueError,
        match="Not valid group_label value",
    ):
        check_group_label(label)


def test_extract_group_ids_empty(tmp_path):
    from clinica.utils.group import extract_group_ids

    assert extract_group_ids(str(tmp_path)) == [""]
    assert extract_group_ids(str(tmp_path / "group")) == [""]
    (tmp_path / "group").mkdir()
    assert extract_group_ids(str(tmp_path / "group")) == [""]


def test_extract_group_ids(tmp_path):
    from clinica.utils.group import extract_group_ids

    groups_folder = tmp_path / "groups"
    groups_folder.mkdir()
    expected_group_ids = ["group-AD", "group-foo", "group-X", "foo"]
    for group_id in expected_group_ids:
        (groups_folder / group_id).mkdir()

    assert set(extract_group_ids(str(tmp_path))) == set(expected_group_ids)

    # The extract_group_ids function currently gets easily fooled
    # as it doesn't perform any checks on the content of the group folder...
    (groups_folder / "file.txt").touch()

    assert set(extract_group_ids(str(tmp_path))) == set(
        expected_group_ids + ["file.txt"]
    )
