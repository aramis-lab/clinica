from dataclasses import dataclass
from os import PathLike
from typing import Union

import pytest


def test_ApplySegmentationDeformation():
    """Uses SPM to apply a deformation field obtained from Segmentation routine to a given file"""

    import clinica.pydra.t1_volume.t1_volume_tasks as t1vol_tasks

    inv = t1vol_tasks.ApplySegmentationDeformation()
    assert inv._jobname == "defs"
    assert inv._jobtype == "util"
    assert inv.input_spec().get_traitsfree() == {
        "mask": 0,
        "mfile": True,
        "use_v8struct": True,
    }
    return
