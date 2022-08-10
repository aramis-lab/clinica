import nest_asyncio
import pydra

nest_asyncio.apply()
from os import getenv
from pathlib import Path

import dwi_preprocessing_using_t1_tasks as dt
import dwi_preprocessing_using_t1_utils as du
import nipype.interfaces.utility as niu
from nipype.interfaces.ants import CreateJacobianDeterminantImage
from nipype.interfaces.ants.registration import RegistrationSynQuick
from nipype.interfaces.c3 import C3dAffineTool
from nipype.interfaces.fsl import MultiImageMaths
from nipype.interfaces.fsl.epi import Eddy
from nipype.interfaces.fsl.maths import Threshold
from nipype.interfaces.fsl.preprocess import BET, FLIRT, ApplyXFM
from nipype.interfaces.fsl.utils import Merge, Split
from nipype.interfaces.io import BIDSDataGrabber
from nipype.interfaces.mrtrix3 import DWIBiasCorrect
from nipype.interfaces.mrtrix3.preprocess import DWIBiasCorrect

bids_data_grabber = BIDSDataGrabber(
    outfields=["t1w", "dwi", "dwi_json", "bvec", "bval"],
    output_query={
        "t1w": {
            "datatype": "anat",
            "suffix": "T1w",
            "extension": ".nii.gz",
        },
        "dwi": {"datatype": "dwi", "suffix": "dwi", "extension": ".nii.gz"},
        "dwi_json": {"datatype": "dwi", "suffix": "dwi", "extension": ".json"},
        "bvec": {"datatype": "dwi", "suffix": "dwi", "extension": ".bvec"},
        "bval": {"datatype": "dwi", "suffix": "dwi", "extension": ".bval"},
    },
)
bet_2 = BET(frac=0.3, robust=True)
use_cuda = True
initrand = False
eddy = Eddy(repol=True, use_cuda=use_cuda, initrand=initrand)
split_epi = Split(dimension="t")
select_epi = niu.Select()
flirt_b0_2_t1 = FLIRT(dof=6, interp="spline", cost="normmi", cost_func="normmi")
ants_reg = RegistrationSynQuick(transform_type="br", dimension=3)
c3d_affine_tool = C3dAffineTool(itk_transform=True, fsl2ras=True)
merge_transform = niu.Merge(3)
jacobian = CreateJacobianDeterminantImage(
    imageDimension=3, outputImage="Jacobian_image.nii.gz"
)
jacmult = MultiImageMaths(op_string="-mul %s")
thres_epi = Threshold(thresh=0.0)
merge_epi = Merge(dimension="t")
bias = DWIBiasCorrect(use_ants=True)
end_bet = BET(mask=True, robust=True)

## Workflow definition
from pydra import Submitter, Workflow
from pydra.tasks.nipype1.utils import Nipype1Task


def build_core_worfflow():
    """Core Workflow  for the DWI preprocessing using T1.

    :param name: The name of the workflow.
    :return: The core workflow.
    """

    workflow = Workflow(
        name="dwi_preprocessing_using_t1",
        # cache_dir = working_dir,
        input_spec=["t1w", "dwi", "dwi_json", "bvec", "bval"],
    )
    # init_input_node
    workflow.add(
        dt.init_input_node(
            name="init_input_node",
            t1w=workflow.lzin.t1w,
            dwi=workflow.lzin.dwi,
            dwi_json=workflow.lzin.dwi_json,
            bvec=workflow.lzin.bvec,
            bval=workflow.lzin.bval,
        )
    )
    workflow.init_input_node.split(("t1w", "dwi", "dwi_json", "bvec", "bval"))
    workflow.add(
        dt.prepare_reference_b0_1(
            name="prepare_reference_b0_1",
            in_dwi=workflow.init_input_node.lzout.dwi,
            in_bval=workflow.init_input_node.lzout.bval,
            in_bvec=workflow.init_input_node.lzout.bvec,
            low_bval=5,
            working_directory=working_dir,
        )
    )
    workflow.add(
        Nipype1Task(
            name="bet_2",
            interface=bet_2,
            in_file=workflow.prepare_reference_b0_1.lzout.out_reference_b0,
            mask=True,
        )
    )
    eddyy = du.build_eddy_wf(
        workflow.prepare_reference_b0_1.lzout.out_b0_dwi_merge,
        workflow.init_input_node.lzout.phase_encoding_direction,
        workflow.init_input_node.lzout.total_readout_time,
        workflow.prepare_reference_b0_1.lzout.out_updated_bval,
        workflow.bet_2.lzout.mask_file,
        workflow.prepare_reference_b0_1.lzout.out_updated_bvec,
    )

    workflow.add(eddyy)

    epii = du.build_epi_wf(
        workflow.init_input_node.lzout.t1w,
        workflow.prepare_reference_b0_1.lzout.out_updated_bvec,
        workflow.init_input_node.lzout.t1w,
    )

    workflow.add(epii)
    workflow.add(
        Nipype1Task(
            name="bias",
            interface=bias,
            in_bval=workflow.prepare_reference_b0_1.lzout.out_updated_bval,
            in_file=workflow.epi.lzout.merge_epi,
            in_bvec=workflow.epi.lzout.rotated_bvecs,
        )
    )
    workflow.add(
        dt.compute_average_b0(
            name="compute_average_b0",
            # interface=compute_average_b0,
            in_dwi=workflow.bias.lzout.out_file,
            in_bval=workflow.prepare_reference_b0_1.lzout.out_updated_bval,
            low_bval=5.0,
        )
    )
    workflow.add(
        Nipype1Task(
            name="end_bet",
            interface=end_bet,
            in_file=workflow.compute_average_b0.lzout.out_b0_average,
        )
    )
    workflow.set_output(
        [
            ("preproc_bval", workflow.prepare_reference_b0_1.lzout.out_updated_bval),
            ("preproc_dwi", workflow.bias.lzout.out_file),
            ("preproc_bvec", workflow.rotate_bvecs.lzout.out_file),
            ("b0_mask", workflow.end_bet.lzout.mask_file),
        ]
    )
    return workflow
