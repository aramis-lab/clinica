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
    from nipype.interfaces.io import BIDSDataGrabber

    from clinica.pydra.interfaces import bids_reader

    task = bids_reader(BIDSQuery(), tmp_path)
    assert isinstance(task, Nipype1Task)
    assert task.name == "bids_reader_task"
    assert task.inputs.base_dir == tmp_path
    assert isinstance(task._interface, BIDSDataGrabber)


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
    "query,grabber",
    [
        (CAPSFileQuery(), CAPSFileDataGrabber),
        (CAPSGroupQuery(), CAPSGroupDataGrabber),
    ],
)
def test_caps_reader_instantiation(tmp_path, query, grabber):
    task = caps_reader(query, tmp_path)
    assert isinstance(task, Nipype1Task)
    assert task.name == "caps_reader_task"
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
        {"mask_tissues": {"tissue_number": (1, 2, 3), "modulation": False}}
    )
    task = caps_reader(query, tmp_path)
    results = run(task)
    assert len(results.output.mask_tissues)
    assert [len(x) for x in results.output.mask_tissues] == [3, 3, 3]
    assert set([Path(x).name for x in results.output.mask_tissues[0]]) == {
        "sub-01_ses-M06_T1w_segm-graymatter_space-Ixi549Space_modulated-off_probability.nii.gz",
        "sub-01_ses-M00_T1w_segm-graymatter_space-Ixi549Space_modulated-off_probability.nii.gz",
        "sub-03_ses-M00_T1w_segm-graymatter_space-Ixi549Space_modulated-off_probability.nii.gz",
    }
    assert set([Path(x).name for x in results.output.mask_tissues[1]]) == {
        "sub-01_ses-M06_T1w_segm-whitematter_space-Ixi549Space_modulated-off_probability.nii.gz",
        "sub-01_ses-M00_T1w_segm-whitematter_space-Ixi549Space_modulated-off_probability.nii.gz",
        "sub-03_ses-M00_T1w_segm-whitematter_space-Ixi549Space_modulated-off_probability.nii.gz",
    }
    assert set([Path(x).name for x in results.output.mask_tissues[2]]) == {
        "sub-01_ses-M06_T1w_segm-csf_space-Ixi549Space_modulated-off_probability.nii.gz",
        "sub-01_ses-M00_T1w_segm-csf_space-Ixi549Space_modulated-off_probability.nii.gz",
        "sub-03_ses-M00_T1w_segm-csf_space-Ixi549Space_modulated-off_probability.nii.gz",
    }


def test_caps_reader_error(tmp_path):
    with pytest.raises(
        TypeError,
        match="caps_reader received an unexpected",
    ):
        caps_reader(BIDSQuery(), tmp_path)  # noqa
