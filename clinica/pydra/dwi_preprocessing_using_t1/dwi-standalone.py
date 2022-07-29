import nest_asyncio
import pydra

nest_asyncio.apply()
from os import getenv
from pathlib import Path, PosixPath

import dwi_task as dt
import toimportfunctions

clinica_in = "/localdrive10TB/users/matthieu.joulot"
pipeline_dir = Path(clinica_in) / "DWIPreprocessingUsingT1"

# pipeline_dir_out = Path(clinica_data_ci_dir) / "DWIPreprocessingUsingT1"
pipeline_in_dir = pipeline_dir / "in"
pipeline_ref_dir = pipeline_dir / "ref"
pipeline_out_dir = pipeline_dir / "out"
working_dir = "/localdrive10TB/users/matthieu.joulot/werk"

# eddy pipeline outputs stored for testing to waste less time
eddy_out_params = "/localdrive10TB/users/matthieu.joulot/DWIPreprocessingUsingT1/eddy/eddy_corrected.eddy_parameters"
eddy_out_corrected = "/localdrive10TB/users/matthieu.joulot/DWIPreprocessingUsingT1/eddy/eddy_corrected.nii.gz"
eddy_out_rotated_bvecs = "/localdrive10TB/users/matthieu.joulot/DWIPreprocessingUsingT1/eddy/eddy_corrected.eddy_rotated_bvecs"

t1w_2 = "/localdrive10TB/users/matthieu.joulot/DWIPreprocessingUsingT1/results1_5/sub-PREVDEMALS0010025PG_ses-M00_T1w.nii.gz"
dwi_2 = "/localdrive10TB/users/matthieu.joulot/DWIPreprocessingUsingT1/results1_5/sub-PREVDEMALS0010025PG_ses-M00_dwi.nii.gz"
dwi_json_2 = "/localdrive10TB/users/matthieu.joulot/DWIPreprocessingUsingT1/results1_5/sub-PREVDEMALS0010025PG_ses-M00_dwi.json"
bval_2 = "/localdrive10TB/users/matthieu.joulot/DWIPreprocessingUsingT1/results1_5/sub-PREVDEMALS0010025PG_ses-M00_dwi.bval"
bvec_2 = "/localdrive10TB/users/matthieu.joulot/DWIPreprocessingUsingT1/results1_5/sub-PREVDEMALS0010025PG_ses-M00_dwi.bvec"
bval_pre = (
    "/localdrive10TB/users/matthieu.joulot/DWIPreprocessingUsingT1/results1_5/bvals"
)
bval_alt = "/localdrive10TB/users/matthieu.joulot/DWIPreprocessingUsingT1/alt_in/bvals"

rotate = "/localdrive10TB/users/matthieu.joulot/rotate/FunctionTask_42137cd3b13d97637d8d68b55680605e6bc6a0e1fda10ef439f2fc97c55b2700/sub-PREVDEMALS0010025PG_ses-M00_dwi_rotated.bvec"

import nipype.interfaces.utility as niu
from nipype.interfaces.ants import CreateJacobianDeterminantImage
from nipype.interfaces.ants.registration import RegistrationSynQuick
from nipype.interfaces.c3 import C3dAffineTool
from nipype.interfaces.fsl import MultiImageMaths
from nipype.interfaces.fsl.epi import Eddy
from nipype.interfaces.fsl.maths import Threshold
from nipype.interfaces.fsl.preprocess import BET, FLIRT
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
use_cuda = False
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


def build_core_worflow():
    """Core Workflow  for the DWI preprocessing using T1."""

    workflow = Workflow(
        name="dwi_preprocessing_using_t1",
        # cache_dir = working_dir,
        input_spec=["input_dir"],
        input_dir=str(pipeline_in_dir / "bids"),
    )

    # bids_data_grabber
    workflow.add(
        Nipype1Task(
            name="bids_data_grabber",
            interface=bids_data_grabber,
            base_dir=workflow.lzin.input_dir,
        )
    )
    # init_input_node
    workflow.add(
        dt.init_input_node(
            name="init_input_node",
            # interface=dt.init_input_node,
            t1w=workflow.bids_data_grabber.lzout.t1w,
            dwi=workflow.bids_data_grabber.lzout.dwi,
            dwi_json=workflow.bids_data_grabber.lzout.dwi_json,
            bvec=workflow.bids_data_grabber.lzout.bvec,
            bval=workflow.bids_data_grabber.lzout.bval,
        ).split(("t1w", "dwi", "dwi_json", "bvec", "bval"))
    )
    workflow.add(
        dt.prepare_reference_b0_1(
            name="prepare_reference_b0_1",
            # interface=dt.prepare_reference_b0_1,
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
    ## Eddy pipeline parts: generate_acq, generate_index, eddy
    workflow.add(
        dt.generate_acq_file(
            name="generate_acq_file",
            # interface=generate_acq_file,
            in_dwi=workflow.prepare_reference_b0_1.lzout.out_b0_dwi_merge,
            fsl_phase_encoding_direction=workflow.init_input_node.lzout.phase_encoding_direction,
            total_readout_time=workflow.init_input_node.lzout.total_readout_time,
            image_id=None,
        )
    )
    workflow.add(
        dt.generate_index_file(
            name="generate_index_file",
            # interface=generate_index_file,
            in_bval=workflow.prepare_reference_b0_1.lzout.out_updated_bval,
            low_bval=5.0,
            image_id=None,
        )
    )
    # workflow.add(
    #     Nipype1Task(
    #         name="eddy",
    #         interface=eddy,
    #         in_file=workflow.prepare_reference_b0_1.lzout.out_b0_dwi_merge,
    #         in_mask=workflow.bet_2.lzout.mask_file,
    #         in_bval=workflow.prepare_reference_b0_1.lzout.out_updated_bval,
    #         in_bvec=workflow.prepare_reference_b0_1.lzout.out_updated_bvec,
    #         in_acqp=workflow.generate_acq_file.lzout.out_acq,
    #         in_index=workflow.generate_index_file.lzout.out_index,
    #     )
    # )
    # eddy pipeline end
    workflow.add(
        Nipype1Task(
            name="split_epi",
            interface=split_epi,
            # in_file=eddy_out_corrected,
            # in_file=workflow.eddy.lzout.out_corrected,
            in_file=workflow.init_input_node.lzout.t1w,
        )
    )
    # workflow2.split("split_epi.lzout.out_files")
    workflow.add(
        Nipype1Task(
            name="select_epi",
            interface=select_epi,
            inlist=workflow.split_epi.lzout.out_files,
            index=[0],
        )
    )
    workflow.add(
        Nipype1Task(
            name="flirt_b0_2_t1",
            interface=flirt_b0_2_t1,
            # in_file=workflow.init_input_node.lzout.t1w,
            in_file=workflow.select_epi.lzout.out,
            reference=workflow.init_input_node.lzout.t1w,
        )
    )
    workflow.add(
        dt.printer(
            name="printer",
            var=workflow.init_input_node.lzout.t1w,
        )
    )
    workflow.add(
        dt.expend_matrix_list(
            name="expend_matrix_list",
            # interface=expend_matrix_list,
            in_matrix=workflow.flirt_b0_2_t1.lzout.out_matrix_file,
            # in_bvec=workflow.eddy.lzout.out_rotated_bvecs,
            in_bvec=workflow.prepare_reference_b0_1.lzout.out_updated_bvec,
        )
    )
    workflow.add(
        dt.rotate_bvecs(
            name="rotate_bvecs",
            # interface=rotate_bvecs,
            # in_bvec=workflow.eddy.lzout.out_rotated_bvecs,
            in_bvec=workflow.prepare_reference_b0_1.lzout.out_updated_bvec,
            in_matrix=workflow.expend_matrix_list.lzout.out_matrix_list,
        )
    )
    workflow.add(
        Nipype1Task(
            name="ants_reg",
            interface=ants_reg,
            moving_image=workflow.flirt_b0_2_t1.lzout.out_file,
            # moving_image=workflow.select_epi.lzout.out,
            # moving_image = workflow.eddy.lzout.out_corrected,
            # fixed_image=workflow.init_input_node.lzout.t1w,
            # fixed_image=workflow.init_input_node.lzout.t1w,
            fixed_image=workflow.printer.lzout.var_out,
        )
    )
    workflow.add(
        Nipype1Task(
            name="c3d_affine_tool",
            interface=c3d_affine_tool,
            source_file=workflow.select_epi.lzout.out,
            # reference_file=workflow.init_input_node.lzout.t1w,
            reference_file=workflow.printer.lzout.var_out,
            transform_file=workflow.flirt_b0_2_t1.lzout.out_matrix_file,
        )
    )
    workflow.add(
        dt.change_itk_transform_type(
            name="change_itk_transform_type",
            # interface=change_itk_transform_type,
            input_affine_file=workflow.c3d_affine_tool.lzout.itk_transform,
        )
    )
    workflow.add(
        Nipype1Task(
            name="merge_transform",
            interface=merge_transform,
            in1=workflow.change_itk_transform_type.lzout.updated_affine_file,
            in2=workflow.ants_reg.lzout.out_matrix,
            in3=workflow.ants_reg.lzout.forward_warp_field,
        )
    )
    workflow.add(
        dt.ants_apply_transforms(
            name="apply_transform_image",
            # interface=ants_apply_transforms,
            fixed_image=workflow.init_input_node.lzout.t1w,
            moving_image=workflow.split_epi.lzout.out_files,
            transforms=workflow.merge_transform.lzout.out,
            warped_image="out_warped_field.nii.gz",
            output_warped_image=True,
        ).split("moving_image")
    )
    workflow.add(
        dt.ants_apply_transforms(
            name="apply_transform_field",
            # interface=ants_apply_transforms,
            fixed_image=workflow.init_input_node.lzout.t1w,
            moving_image=workflow.split_epi.lzout.out_files,
            transforms=workflow.merge_transform.lzout.out,
            warped_image="out_warped.nii.gz",
            output_warped_image=False,
        ).split("moving_image")
    )
    workflow.add(
        Nipype1Task(
            name="jacobian",
            interface=jacobian,
            deformationField=workflow.apply_transform_field.lzout.warped_image,
        )
    )
    workflow.jacobian.combine(["apply_transform_field.moving_image"])
    workflow.apply_transform_image.combine(["moving_image"])
    workflow.add(
        Nipype1Task(
            name="jacmult",
            interface=jacmult,
            operand_files=workflow.apply_transform_image.lzout.warped_image,
            in_file=workflow.jacobian.lzout.jacobian_image,
        ).split(("operand_files", "in_file"))
    )
    workflow.add(
        Nipype1Task(
            name="thres_epi",
            interface=thres_epi,
            in_file=workflow.jacmult.lzout.out_file,
        )
    )
    workflow.thres_epi.combine(["jacmult.operand_files", "jacmult.in_file"])

    workflow.add(
        Nipype1Task(
            name="merge_epi",
            interface=merge_epi,
            in_files=workflow.thres_epi.lzout.out_file,
        )
    )
    ## epi pipeline end
    workflow.add(
        Nipype1Task(
            name="bias",
            interface=bias,
            in_bval=workflow.prepare_reference_b0_1.lzout.out_updated_bval,
            in_file=workflow.merge_epi.lzout.merged_file,
            in_bvec=workflow.rotate_bvecs.lzout.out_file,
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
            # ("bdg", workflow.bids_data_grabber.lzout.t1w),
            # ("init", workflow.init_input_node.lzout.t1w,),
            # ("dvfv", workflow.split_epi.lzout.out_files),
            # ("select", workflow.select_epi.lzout.out),
            # ("ants_registration_1", workflow.ants_reg.lzout.forward_warp_field),
            # ("eddy_out_corrected", workflow.eddy.lzout.out_corrected),
            # ("flirt", workflow.flirt_b0_2_t1.lzout.out_file),
            # ("t1w", workflow.init_input_node.lzout.t1w),
            # ("preproc_bval", workflow.prepare_reference_b0_1.lzout.out_updated_bval),
            # ("preproc_dwi", workflow.bias.lzout.out_file),
            # ("preproc_bvec", workflow.rotate_bvecs.lzout.out_file),
            # ("b0_mask", workflow.end_bet.lzout.mask_file),
            # ("c3d", workflow.c3d_affine_tool.lzout.itk_transform),
            # ("expand", workflow.expend_matrix_list.lzout.out_matrix_list),
            # final results to use for after
            ("preproc_dwi", workflow.bias.lzout.out_file),
            ("preproc_bvec", workflow.rotate_bvecs.lzout.out_file),
            ("preproc_bval", workflow.prepare_reference_b0_1.lzout.out_updated_bval),
            ("b0_mask", workflow.end_bet.lzout.mask_file),
            # ("mask_file",workflow.bet_2.lzout.mask_file),
            # ("in_file",workflow.prepare_reference_b0_1.lzout.out_b0_dwi_merge),
            # ("in_bval",workflow.prepare_reference_b0_1.lzout.out_updated_bval),
            # ("in_bvec",workflow.prepare_reference_b0_1.lzout.out_updated_bvec),
            # ("out_acq",workflow.generate_acq_file.lzout.out_acq),
            # ("out_index",workflow.generate_index_file.lzout.out_index),
        ]
    )
    return workflow
    # in_file=workflow.prepare_reference_b0_1.lzout.out_b0_dwi_merge,
    #         in_mask=workflow.bet_2.lzout.mask_file,
    #         in_bval=workflow.prepare_reference_b0_1.lzout.out_updated_bval,
    #         in_bvec=workflow.prepare_reference_b0_1.lzout.out_updated_bvec,
    #         in_acqp=workflow.generate_acq_file.lzout.out_acq,
    #         in_index=workflow.generate_index_file.lzout.out_index,


def run_wf(workflow):
    with Submitter(plugin="cf") as submitter:
        submitter(workflow)
    results = workflow.result(return_inputs=True)
    print(results)
    return results


def build_run():
    workflow = build_core_worflow()
    results = run_wf(workflow)

    return results


results = build_run()

# print("\nresults:\n", results)
print("\n\n\n\nresults2:\n", results)
