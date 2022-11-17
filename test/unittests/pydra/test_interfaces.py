from pathlib import Path

import pytest
from pydra.tasks.nipype1.utils import Nipype1Task

from clinica.pydra.interfaces import (
    CAPSFileDataGrabber,
    CAPSGroupDataGrabber,
    caps_reader,
)
from clinica.pydra.query import BIDSQuery, CAPSFileQuery, CAPSGroupQuery
from clinica.utils.testing_utils import build_bids_directory, build_caps_directory


def test_bids_reader_instantiation(tmp_path):
    from clinica.pydra.interfaces import bids_reader
    from pydra.engine.task import FunctionTask

    task = bids_reader(BIDSQuery(), tmp_path)
    assert isinstance(task, FunctionTask)
    assert task.name == "bids_reader_task"
    assert task.inputs.dataset_path == tmp_path


def test_bids_reader(tmp_path):
    from clinica.pydra.engine_utils import run
    from clinica.pydra.interfaces import bids_reader

    structure = {
        "sub-01": ["ses-M00", "ses-M06"],
        "sub-03": ["ses-M00"],
    }
    build_bids_directory(tmp_path, structure)
    query = BIDSQuery({"T1w": {}})
    task = bids_reader(query, tmp_path)
    results = run(task)
    assert set([Path(f).name for f in results.output.T1w]) == {
        "sub-01_ses-M00_T1w.nii.gz",
        "sub-01_ses-M06_T1w.nii.gz",
        "sub-03_ses-M00_T1w.nii.gz",
    }


@pytest.mark.parametrize(
    "query,grabber,name",
    [
        (CAPSFileQuery(), CAPSFileDataGrabber, "caps_file_reader_task"),
        (CAPSGroupQuery(), CAPSGroupDataGrabber, "caps_group_reader_task"),
    ],
)
def test_caps_reader_instantiation(tmp_path, query, grabber, name):
    task = caps_reader(query, tmp_path)
    assert isinstance(task, Nipype1Task)
    assert task.name == name
    assert task.inputs.base_dir == tmp_path
    assert isinstance(task._interface, grabber)


def test_caps_reader(tmp_path):
    from clinica.pydra.engine_utils import run
    from clinica.pydra.interfaces import caps_reader

    structure = {
        "groups": ["UnitTest"],
        "pipelines": ["t1"],
        "subjects": {
            "sub-01": ["ses-M00", "ses-M06"],
            "sub-03": ["ses-M00"],
        },
    }
    build_caps_directory(tmp_path, structure)
    query = CAPSFileQuery(
        {
            "mask_tissues": {"tissue_number": (1, 2, 3), "modulation": False},
            "flow_fields": {"group_label": "UnitTest"},
            "pvc_mask_tissues": {"tissue_number": (1, 2, 3)},
            "dartel_template": {"group_label": "UnitTest"},
        }
    )
    task = caps_reader(query, tmp_path)
    results = run(task)
    for key, _ in query.query.items():
        assert hasattr(results.output, key)
    for key in ["mask_tissues", "pvc_mask_tissues"]:
        assert len(getattr(results.output, key)) == 3
        assert [len(x) for x in getattr(results.output, key)] == [3, 3, 3]
        space_tag = "_space-Ixi549Space_modulated-off" if key == "mask_tissues" else ""
        for i, tissue in enumerate(["graymatter", "whitematter", "csf"]):
            assert set([Path(x).name for x in getattr(results.output, key)[i]]) == {
                f"sub-01_ses-M06_T1w_segm-{tissue}{space_tag}_probability.nii.gz",
                f"sub-01_ses-M00_T1w_segm-{tissue}{space_tag}_probability.nii.gz",
                f"sub-03_ses-M00_T1w_segm-{tissue}{space_tag}_probability.nii.gz",
            }
    assert len(results.output.flow_fields) == 3
    assert set([Path(x).name for x in results.output.flow_fields]) == {
        "sub-01_ses-M06_T1w_target-UnitTest_transformation-forward_deformation.nii.gz",
        "sub-01_ses-M00_T1w_target-UnitTest_transformation-forward_deformation.nii.gz",
        "sub-03_ses-M00_T1w_target-UnitTest_transformation-forward_deformation.nii.gz",
    }
    query = CAPSGroupQuery({"dartel_template": {"group_label": "UnitTest"}})
    task = caps_reader(query, tmp_path)
    results = run(task)
    assert hasattr(results.output, "dartel_template")
    assert Path(results.output.dartel_template).name == "group-UnitTest_template.nii.gz"


def test_caps_reader_error(tmp_path):
    with pytest.raises(
        TypeError,
        match="caps_reader received an unexpected",
    ):
        caps_reader(BIDSQuery(), tmp_path)  # noqa
