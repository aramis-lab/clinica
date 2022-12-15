import pytest


@pytest.mark.parametrize(
    "flowfields,tissues,expected",
    [
        (
            ["path_to_flowfield_1", "path_to_flowfield_2"],
            (1, 2, 3),
            [['path_to_flowfield_1'] * 3, ['path_to_flowfield_2'] * 3],
        ),
        (
            ["path_to_flowfield_1", "path_to_flowfield_2"],
            (1, 3),
            [['path_to_flowfield_1'] * 2, ['path_to_flowfield_2'] * 2],
        ),
        (
            ["path_to_flowfield_1"],
            (1, 3),
            [['path_to_flowfield_1'] * 2],
        ),
        (
            "path_to_flowfield_1",
            (1, 3),
            ['path_to_flowfield_1'] * 2,
        ),
        (
            "path_to_flowfield_1",
            "abc",  # This does not make sense but will not raise
            ['path_to_flowfield_1'] * 3,
        ),
    ]
)
def test_prepare_flowfields(flowfields, tissues, expected):
    from clinica.pydra.t1_volume.dartel2mni.tasks import prepare_flowfields

    assert prepare_flowfields(flowfields, tissues) == expected


def test_prepare_flowfields_errors():
    from clinica.pydra.t1_volume.dartel2mni.tasks import prepare_flowfields
    
    with pytest.raises(
        ValueError,
        match="Unvalid tissues type: <class 'int'>.",
    ):
        prepare_flowfields("path", 1)
    
    with pytest.raises(
        ValueError,
        match="Unvalid flowfields type: <class 'int'>.",
    ):
        prepare_flowfields(42, (1, 2))
