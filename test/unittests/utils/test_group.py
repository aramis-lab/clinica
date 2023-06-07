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


@pytest.mark.parametrize("folder_name", ["group-AD", "group-foo", "group-foo12bar"])
def test_check_group_dir(tmp_path, folder_name):
    from clinica.utils.group import _check_group_dir

    assert _check_group_dir(tmp_path / folder_name) is None


@pytest.mark.parametrize(
    "folder_name", ["foo", "foo-group", "groupfoo12bar", "foo.txt"]
)
def test_check_group_dir_errors(tmp_path, folder_name):
    from clinica.utils.group import _check_group_dir

    with pytest.raises(
        ValueError,
        match="Group directory",
    ):
        _check_group_dir(tmp_path / folder_name)


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
        ValueError, match=f"Group directory {groups_folder / 'foo'} is not valid"
    ):
        extract_group_ids(tmp_path)
