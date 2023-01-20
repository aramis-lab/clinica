import pytest


@pytest.mark.parametrize(
    "flowfields,number_of_tissues,expected",
    [
        (
            ["path_to_flowfield_1", "path_to_flowfield_2"],
            3,
            [["path_to_flowfield_1"] * 3, ["path_to_flowfield_2"] * 3],
        ),
        (
            ["path_to_flowfield_1", "path_to_flowfield_2"],
            2,
            [["path_to_flowfield_1"] * 2, ["path_to_flowfield_2"] * 2],
        ),
        (
            ["path_to_flowfield_1"],
            2,
            [["path_to_flowfield_1"] * 2],
        ),
        (
            "path_to_flowfield_1",
            2,
            ["path_to_flowfield_1"] * 2,
        ),
    ],
)
def test_prepare_flowfields(flowfields, number_of_tissues, expected):
    from clinica.pydra.t1_volume.dartel2mni.tasks import prepare_flowfields

    assert prepare_flowfields(flowfields, number_of_tissues) == expected


def test_prepare_flowfields_errors():
    from clinica.pydra.t1_volume.dartel2mni.tasks import prepare_flowfields

    with pytest.raises(
        TypeError,
        match="Invalid type for number_of_tissues. Expected int, got <class 'str'>.",
    ):
        prepare_flowfields("path", "abc")

    with pytest.raises(
        TypeError,
        match="Invalid flowfields type: <class 'int'>.",
    ):
        prepare_flowfields(42, 2)
