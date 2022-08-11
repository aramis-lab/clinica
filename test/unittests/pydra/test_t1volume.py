from pydra import Workflow
import pytest
from os import PathLike

from typing import Union


@pytest.fixture
def nifti_list():
    path_in = "/Users/omar.elrifai/workspace/data_ci/CenterNifti/in/bids/\
        sub-11/ses-M0/anat/sub-11_ses-M0_T1w.nii"
    return path_in


# def test_zip_nii(in_var: Union[PathLike, list], same_dir: bool = False) -> PathLike:
#     """Zips a file or a list of files"""

#     from clinica.pydra.t1_volume.t1_volume_tasks import zip_nii

#     out_file = zip_nii(in_var, False)
#     return out_file


def test_check_volume_location_in_world_coordinate_system(nifti_list):

    from clinica.pydra.t1_volume.t1_volume_tasks import (
        check_volume_location_in_world_coordinate_system,
    )

    task_check_center = check_volume_location_in_world_coordinate_system(
        nifti_list=nifti_list,
        bids_dir="/tmp",
    )

    res = task_check_center()
    assert res.output.out_status == False

    return


# def test_ApplySegmentationDeformation():
#     """Uses SPM to apply a deformation field obtained from Segmentation routine to a given file"""

#     import clinica.pydra.t1_volume.t1_volume_tasks as t1vol_tasks

#     inv = t1vol_tasks.ApplySegmentationDeformation()
#     assert inv._jobname == "defs"
#     assert inv._jobtype == "util"
#     assert inv.input_spec().get_traitsfree() == {
#         "mask": 0,
#         "mfile": True,
#         "use_v8struct": True,
#     }
#     return

if __name__ == "__main__":
    test_check_volume_location_in_world_coordinate_system(nifti_list)
