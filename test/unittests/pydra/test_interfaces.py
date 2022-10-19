import pytest
from pydra.tasks.nipype1.utils import Nipype1Task

from clinica.pydra.interfaces import (
    CAPSFileDataGrabber,
    CAPSGroupDataGrabber,
    caps_reader,
)
from clinica.pydra.query import BIDSQuery, CAPSFileQuery, CAPSGroupQuery


def test_bids_reader(tmp_path):
    from nipype.interfaces.io import BIDSDataGrabber

    from clinica.pydra.interfaces import bids_reader

    task = bids_reader(BIDSQuery(), tmp_path)
    assert isinstance(task, Nipype1Task)
    assert task.name == "bids_reader_task"
    assert task.inputs.base_dir == tmp_path
    assert isinstance(task._interface, BIDSDataGrabber)


@pytest.mark.parametrize(
    "query,grabber",
    [
        (CAPSFileQuery(), CAPSFileDataGrabber),
        (CAPSGroupQuery(), CAPSGroupDataGrabber),
    ],
)
def test_caps_reader(tmp_path, query, grabber):
    task = caps_reader(query, tmp_path)
    assert isinstance(task, Nipype1Task)
    assert task.name == "caps_reader_task"
    assert task.inputs.base_dir == tmp_path
    assert isinstance(task._interface, grabber)


def test_caps_reader_error(tmp_path):
    with pytest.raises(
        TypeError,
        match="caps_reader received an unexpected",
    ):
        caps_reader(BIDSQuery(), tmp_path)  # noqa
