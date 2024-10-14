from pathlib import Path

import pytest

from clinica.pydra.query import BIDSQuery, CAPSFileQuery, CAPSGroupQuery, Query
from clinica.utils.input_files import QueryPattern


def test_query():
    q = Query()
    assert q.query == {}
    assert len(q) == 0


def test_bids_query():
    q = BIDSQuery()
    assert q.query == {}
    assert len(q) == 0

    q = BIDSQuery({"T1w": {}})
    assert q.query == {"T1w": {"suffix": "T1w", "extension": [".nii.gz"]}}
    assert len(q) == 1

    q = BIDSQuery({"T1w": {"foo": "bar"}})
    assert q.query == {
        "T1w": {
            "suffix": "T1w",
            "extension": [".nii.gz"],
            "foo": "bar",
        }
    }
    assert len(q) == 1

    q = BIDSQuery({"T1w": {}, "foo": {"bar": "baz"}})
    assert q.query == {"T1w": {"suffix": "T1w", "extension": [".nii.gz"]}}
    assert len(q) == 1


def test_caps_file_query():
    q = CAPSFileQuery({"mask_tissues": {"tissue_number": (1, 2), "modulation": False}})
    assert len(q) == 1
    assert q.query == {
        "mask_tissues": [
            QueryPattern(
                str(
                    Path("t1")
                    / "spm"
                    / "segmentation"
                    / "normalized_space"
                    / "*_*_T1w_segm-graymatter_space-Ixi549Space_modulated-off_probability.nii*"
                ),
                "Tissue probability map graymatter based on native MRI in MNI space (Ixi549) without modulation.",
                "t1-volume-tissue-segmentation",
            ),
            QueryPattern(
                str(
                    Path("t1")
                    / "spm"
                    / "segmentation"
                    / "normalized_space"
                    / "*_*_T1w_segm-whitematter_space-Ixi549Space_modulated-off_probability.nii*"
                ),
                "Tissue probability map whitematter based on native MRI in MNI space (Ixi549) without modulation.",
                "t1-volume-tissue-segmentation",
            ),
        ]
    }

    q = CAPSFileQuery(
        {
            "mask_tissues": {"tissue_number": (1, 2), "modulation": False},
            "flow_fields": {"group_label": "UnitTest"},
        }
    )
    assert len(q) == 2
    assert q.query == {
        "mask_tissues": [
            QueryPattern(
                str(
                    Path("t1")
                    / "spm"
                    / "segmentation"
                    / "normalized_space"
                    / "*_*_T1w_segm-graymatter_space-Ixi549Space_modulated-off_probability.nii*"
                ),
                "Tissue probability map graymatter based on native MRI in MNI space (Ixi549) without modulation.",
                "t1-volume-tissue-segmentation",
            ),
            QueryPattern(
                str(
                    Path("t1")
                    / "spm"
                    / "segmentation"
                    / "normalized_space"
                    / "*_*_T1w_segm-whitematter_space-Ixi549Space_modulated-off_probability.nii*"
                ),
                "Tissue probability map whitematter based on native MRI in MNI space (Ixi549) without modulation.",
                "t1-volume-tissue-segmentation",
            ),
        ],
        "flow_fields": QueryPattern(
            str(
                Path("t1")
                / "spm"
                / "dartel"
                / "group-UnitTest"
                / "sub-*_ses-*_T1w_target-UnitTest_transformation-forward_deformation.nii*"
            ),
            "Deformation from native space to group template UnitTest space.",
            "t1-volume-create-dartel",
        ),
    }


def test_caps_group_query():
    q = CAPSGroupQuery({"dartel_template": {"group_label": "UnitTest"}})
    assert len(q) == 1
    assert q.query == {
        "dartel_template": QueryPattern(
            str(Path("group-UnitTest") / "t1" / "group-UnitTest_template.nii*"),
            "T1w template file of group UnitTest",
            "t1-volume or t1-volume-create-dartel",
        )
    }


@pytest.mark.parametrize(
    "query_type,expected",
    [("bids", BIDSQuery), ("caps_file", CAPSFileQuery), ("caps_group", CAPSGroupQuery)],
)
def test_query_factory(query_type, expected):
    from clinica.pydra.query import query_factory

    assert isinstance(query_factory({}, query_type), expected)


def test_query_factory_error():
    from clinica.pydra.query import query_factory

    with pytest.raises(ValueError):
        query_factory({}, "fooo")
