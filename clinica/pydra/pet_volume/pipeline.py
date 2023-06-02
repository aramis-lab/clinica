from typing import Any

import pydra
from nipype.algorithms.misc import Gunzip
from nipype.interfaces.petpvc import PETPVC
from nipype.interfaces.spm import Coregister, DARTELNorm2MNI, Smooth
from nipype.interfaces.spm.utils import Reslice
from pydra.engine import Workflow
from pydra.tasks.nipype1.utils import Nipype1Task

from clinica.pydra.engine import clinica_io
from clinica.pydra.pet_volume.tasks import (
    apply_binary_mask_task,
    atlas_statistics_task,
    create_binary_mask_task,
    create_pvc_mask_task,
    get_psf_task,
    normalize_to_reference_task,
    pet_pvc_name_task,
)


def _check_pipeline_parameters(parameters: dict) -> dict:
    """Check the parameters passed to the pipeline.

    Parameters
    ----------
    parameters : dict
        Dictionary of parameters to analyze.

    Returns
    -------
    dict :
        Cleaned dictionary of parameters.
    """
    from clinica.utils.atlas import PET_VOLUME_ATLASES
    from clinica.utils.group import check_group_label

    parameters.setdefault("group_label", None)
    check_group_label(parameters["group_label"])
    if "acq_label" not in parameters.keys():
        raise KeyError("Missing compulsory acq_label key in pipeline parameter.")
    parameters.setdefault("pvc_psf_tsv", None)
    parameters["apply_pvc"] = parameters["pvc_psf_tsv"] is not None
    parameters.setdefault("mask_tissues", [1, 2, 3])
    parameters.setdefault("mask_threshold", 0.3)
    parameters.setdefault("pvc_mask_tissues", [1, 2, 3])
    parameters.setdefault("smooth", 8.0)
    parameters.setdefault("atlases", PET_VOLUME_ATLASES)
    return parameters


@clinica_io
def build_core_workflow(name: str = "core", parameters: dict = {}) -> Workflow:
    """Build the core workflow for the PET Volume pipeline.

    Parameters
    ----------
    name : str, optional
        The name of the workflow. Default="core".

    parameters : dict, optional
        Optional dictionary of parameters to be used
        within the workflow.
        Default={}.

    Returns
    -------
    wf : Workflow
        The core workflow.
    """
    from pydra.tasks import petpvc

    from clinica.pydra.shared_workflows.smoothing import build_smoothing_workflow
    from clinica.pydra.utils import sanitize_fwhm
    from clinica.utils.pet import get_suvr_mask
    from clinica.utils.spm import spm_standalone_is_available, use_spm_standalone

    if spm_standalone_is_available():
        use_spm_standalone()

    parameters = _check_pipeline_parameters(parameters)

    input_spec = pydra.specs.SpecInfo(
        name="Input",
        fields=[
            ("_graph_checksums", Any),
            ("T1w", dict, {}, {"mandatory": True}),
            ("pet", dict, {"tracer": parameters["acq_label"]}, {"mandatory": True}),
            (
                "mask_tissues",
                dict,
                {"tissue_number": parameters["mask_tissues"], "modulation": False},
                {"mandatory": True},
            ),
            (
                "flow_fields",
                dict,
                {"group_label": parameters["group_label"]},
                {"mandatory": True},
            ),
            (
                "dartel_template",
                dict,
                {"group_label": parameters["group_label"]},
                {"mandatory": True},
            ),
            (
                "pvc_mask_tissues",
                dict,
                {"tissue_number": parameters["pvc_mask_tissues"]},
                {"mandatory": False},
            ),
        ],
        bases=(pydra.specs.BaseSpec,),
    )

    wf = Workflow(
        name,
        input_spec=input_spec,
    )
    wf.split(("pet", "T1w", "flow_fields", "pvc_mask_tissues", "mask_tissues"))
    compressed_inputs = [
        "pet",
        "T1w",
        "flow_fields",
        "dartel_template",
    ]
    # Unzipping
    for input_name in compressed_inputs:
        wf.add(
            Nipype1Task(
                name=f"unzip_{input_name}",
                interface=Gunzip(),
                in_file=getattr(wf.lzin, input_name),
            )
        )
    wf.add(
        Nipype1Task(
            name="unzip_mask_tissues",
            interface=Gunzip(),
            in_file=wf.lzin.mask_tissues,
        )
        .split("in_file")
        .combine("in_file")
    )
    wf.add(
        Nipype1Task(
            name="unzip_reference_mask",
            interface=Gunzip(),
            in_file=get_suvr_mask(parameters["suvr_reference_region"]),
        )
    )

    # Coregister PET into T1 native space
    coreg_pet_t1 = Nipype1Task(
        name="coreg_pet_t1",
        interface=Coregister(),
    )
    coreg_pet_t1.inputs.source = wf.unzip_pet.lzout.out_file
    coreg_pet_t1.inputs.target = wf.unzip_T1w.lzout.out_file
    wf.add(coreg_pet_t1)

    # Spatially normalize PET into MNI
    dartel_mni_reg = Nipype1Task(
        name="dartel_mni_reg",
        interface=DARTELNorm2MNI(),
    )
    dartel_mni_reg.inputs.modulate = False
    dartel_mni_reg.inputs.fwhm = 0
    dartel_mni_reg.inputs.flowfield_files = wf.unzip_flow_fields.lzout.out_file
    dartel_mni_reg.inputs.template_file = wf.unzip_dartel_template.lzout.out_file
    dartel_mni_reg.inputs.apply_to_files = wf.coreg_pet_t1.lzout.coregistered_source
    wf.add(dartel_mni_reg)

    # Reslice reference region mask into PET
    reslice = Nipype1Task(
        name="reslice",
        interface=Reslice(),
    )
    reslice.inputs.in_file = wf.unzip_reference_mask.lzout.out_file
    reslice.inputs.space_defining = wf.dartel_mni_reg.lzout.normalized_files
    wf.add(reslice)

    # Normalize PET values according to reference region
    wf.add(
        normalize_to_reference_task(
            name="norm_to_ref",
            interface=normalize_to_reference_task,
            pet_image=wf.dartel_mni_reg.lzout.normalized_files,
            region_mask=wf.reslice.lzout.out_file,
        )
    )

    # Create binary mask from segmented tissues
    wf.add(
        create_binary_mask_task(
            name="binary_mask",
            interface=create_binary_mask_task,
            tissues=wf.unzip_mask_tissues.lzout.out_file,
            threshold=parameters["mask_threshold"],
        )
    )

    # Mask PET image
    wf.add(
        apply_binary_mask_task(
            name="apply_mask",
            interface=apply_binary_mask_task,
            image=wf.norm_to_ref.lzout.suvr_pet_path,
            binary_mask=wf.binary_mask.lzout.out_mask,
        )
    )

    # Smoothing
    if parameters["smooth"] is not None:
        smoothing_wf = build_smoothing_workflow(
            "smoothing_workflow",
            sanitize_fwhm(parameters["smooth"]),
        )
        smoothing_wf.inputs.input_file = wf.apply_mask.lzout.masked_image_path
        wf.add(smoothing_wf)

    # Atlas Statistics
    wf.add(
        atlas_statistics_task(
            name="atlas_stats_node",
            interface=atlas_statistics_task,
            in_image=wf.norm_to_ref.lzout.suvr_pet_path,
            in_atlas_list=parameters["atlases"],
        )
    )

    # PVC
    if parameters["apply_pvc"]:

        # Unzipping
        wf.add(
            Nipype1Task(
                name="unzip_pvc_mask_tissues",
                interface=Gunzip(),
                in_file=wf.lzin.pvc_mask_tissues,
            )
            .split("in_file")
            .combine("in_file")
        )

        # Creating Mask to use in PVC
        wf.add(
            create_pvc_mask_task(
                name="pvc_mask",
                interface=create_pvc_mask_task,
                tissues=wf.unzip_pvc_mask_tissues.lzout.out_file,
            )
        )

        # Build PET PVC name for PET PVC task
        wf.add(
            pet_pvc_name_task(
                name="pet_pvc_name",
                interface=pet_pvc_name_task,
                pet_image=wf.coreg_pet_t1.lzout.coregistered_source,
                pvc_method="RBV",
            )
        )

        # Retrieve PSF for PET PVC task
        wf.add(
            get_psf_task(
                name="get_psf_task",
                interface=get_psf_task,
                pvc_psf_tsv=parameters["pvc_psf_tsv"],
                pet_filename=wf.lzin.pet,
                pet_tracer=parameters["acq_label"],
            )
        )

        # PET PVC
        wf.add(
            petpvc.PETPVC(
                name="pvc",
                pvc_method="RBV",
                output_image=wf.pet_pvc_name.lzout.pet_pvc_path,
                input_image=wf.coreg_pet_t1.lzout.coregistered_source,
                mask_image=wf.pvc_mask.lzout.out_mask,
                fwhm_x=wf.get_psf_task.lzout.psf_x,
                fwhm_y=wf.get_psf_task.lzout.psf_y,
                fwhm_z=wf.get_psf_task.lzout.psf_z,
            )
        )

        # Spatially normalize PET into MNI
        dartel_mni_reg_pvc = Nipype1Task(
            name="dartel_mni_reg_pvc",
            interface=DARTELNorm2MNI(),
        )
        dartel_mni_reg_pvc.inputs.modulate = False
        dartel_mni_reg_pvc.inputs.fwhm = 0
        dartel_mni_reg_pvc.inputs.flowfield_files = wf.unzip_flow_fields.lzout.out_file
        dartel_mni_reg_pvc.inputs.template_file = (
            wf.unzip_dartel_template.lzout.out_file
        )
        dartel_mni_reg_pvc.inputs.apply_to_files = wf.pvc.lzout.output_image
        wf.add(dartel_mni_reg_pvc)

        # Reslice reference region mask into PET
        reslice_pvc = Nipype1Task(
            name="reslice_pvc",
            interface=Reslice(),
        )
        reslice_pvc.inputs.space_defining = wf.dartel_mni_reg_pvc.lzout.normalized_files
        reslice_pvc.inputs.in_file = wf.unzip_reference_mask.lzout.out_file
        wf.add(reslice_pvc)

        # Normalize PET values according to reference region
        wf.add(
            normalize_to_reference_task(
                name="norm_to_ref_pvc",
                interface=normalize_to_reference_task,
                pet_image=wf.dartel_mni_reg_pvc.lzout.normalized_files,
                region_mask=wf.reslice_pvc.lzout.out_file,
            )
        )

        # Mask PET image
        wf.add(
            apply_binary_mask_task(
                name="apply_mask_pvc",
                interface=apply_binary_mask_task,
                image=wf.norm_to_ref_pvc.lzout.suvr_pet_path,
                binary_mask=wf.binary_mask.lzout.out_mask,
            )
        )

        # Smoothing
        if parameters["smooth"] is not None:
            pvc_smoothing_wf = build_smoothing_workflow(
                "pvc_smoothing_workflow",
                sanitize_fwhm(parameters["smooth"]),
            )
            pvc_smoothing_wf.inputs.input_file = (
                wf.apply_mask_pvc.lzout.masked_image_path
            )
            wf.add(pvc_smoothing_wf)

        # Atlas Statistics
        wf.add(
            atlas_statistics_task(
                name="atlas_stats_pvc",
                interface=atlas_statistics_task,
                in_image=wf.norm_to_ref_pvc.lzout.suvr_pet_path,
                in_atlas_list=parameters["atlases"],
            )
        )

    output_connections = [
        ("pet_suvr_masked_smoothed", wf.smoothing_workflow.lzout.smoothed_files),
        ("pet_t1_native", wf.coreg_pet_t1.lzout.coregistered_source),
        ("pet_mni", wf.dartel_mni_reg.lzout.normalized_files),
        ("pet_suvr", wf.norm_to_ref.lzout.suvr_pet_path),
        ("binary_mask", wf.binary_mask.lzout.out_mask),
        ("pet_suvr_masked", wf.apply_mask.lzout.masked_image_path),
        ("atlas_statistics", wf.atlas_stats_node.lzout.atlas_statistics),
    ]

    if parameters["apply_pvc"]:
        output_connections += [
            (
                "pet_pvc_suvr_masked_smoothed",
                wf.pvc_smoothing_workflow.lzout.smoothed_files,
            ),
            ("pet_pvc", wf.pvc.lzout.output_image),
            ("pet_pvc_mni", wf.dartel_mni_reg_pvc.lzout.normalized_files),
            ("pet_pvc_suvr", wf.norm_to_ref_pvc.lzout.suvr_pet_path),
            ("pet_pvc_suvr_masked", wf.apply_mask_pvc.lzout.masked_image_path),
            ("pvc_atlas_statistics", wf.atlas_stats_pvc.lzout.atlas_statistics),
        ]

    wf.set_output(output_connections)

    return wf
