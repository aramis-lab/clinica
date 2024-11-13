import re

import pytest


@pytest.mark.parametrize(
    "label",
    ["A", "a1", "123", "0xY2", "mygroup"],
)
def test_check_group_label(label):
    from clinica.utils.group import GroupLabel

    label = GroupLabel(label)

    assert label == label


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
    from clinica.utils.group import GroupLabel

    with pytest.raises(
        ValueError,
        match=re.escape(
            f"Group label '{label}' is not a valid group label: it must be composed only by letters and/or numbers."
        ),
    ):
        GroupLabel(label)


@pytest.mark.parametrize("folder_name", ["group-AD", "group-foo", "group-foo12bar"])
def test_group_id(folder_name: str):
    from clinica.utils.group import GroupID

    assert GroupID(folder_name) == folder_name


@pytest.mark.parametrize(
    "folder_name", ["foo bar", "", "foo$", "foo-group", "group_foo_12_bar", "foo.txt"]
)
def test_group_id_errors(folder_name):
    from clinica.utils.group import GroupID

    with pytest.raises(
        ValueError,
        match=re.escape(
            f"Group ID '{folder_name}' is not a valid group ID: it must start with 'group-'."
        ),
    ):
        GroupID(folder_name)


def test_extract_group_ids_empty(tmp_path):
    from clinica.utils.group import extract_group_ids

    assert extract_group_ids(tmp_path) == []
    assert extract_group_ids(tmp_path / "group") == []
    (tmp_path / "group").mkdir()
    assert extract_group_ids(tmp_path / "group") == []


def test_extract_group_ids(tmp_path):
    from clinica.utils.group import extract_group_ids

    groups_folder = tmp_path / "groups"
    groups_folder.mkdir()
    expected_group_ids = ["group-AD", "group-foo", "group-X"]
    for group_id in expected_group_ids:
        (groups_folder / group_id).mkdir()

    assert set(extract_group_ids(tmp_path)) == set(expected_group_ids)

    (groups_folder / "file.txt").touch()

    assert set(extract_group_ids(tmp_path)) == set(expected_group_ids)


def test_extract_group_ids_errors(tmp_path):
    from clinica.utils.group import extract_group_ids

    groups_folder = tmp_path / "groups"
    groups_folder.mkdir()
    expected_group_ids = ["group-AD", "group-foo", "group-X", "foo"]
    for group_id in expected_group_ids:
        (groups_folder / group_id).mkdir()

    with pytest.raises(
        ValueError,
        match="Group ID 'foo' is not a valid group ID: it must start with 'group-'.",
    ):
        extract_group_ids(tmp_path)
