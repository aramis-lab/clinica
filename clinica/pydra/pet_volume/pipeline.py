from nipype.interfaces.petpvc import PETPVC
from nipype.interfaces.spm import Coregister, DARTELNorm2MNI, Smooth
from nipype.interfaces.spm.utils import Reslice
from pydra import Workflow
from pydra.tasks.nipype1.utils import Nipype1Task

from clinica.pydra.engine import clinica_io
from clinica.pydra.pet_volume.tasks import (
    apply_binary_mask,
    atlas_statistics,
    create_binary_mask,
    create_pvc_mask,
    normalize_to_reference,
)


@clinica_io
def build_core_workflow(name: str = "core") -> Workflow:
    """Core workflow for the PET Volume pipeline.

    Parameters
    ----------
    name : The name of the workflow.

    Returns
    -------
    workflow : The core workflow.
    """
    from clinica.utils.spm import spm_standalone_is_available, use_spm_standalone

    if spm_standalone_is_available():
        use_spm_standalone()

    wf = Workflow(name, input_spec=["T1w", "pet"])
    compressed_inputs = [
        "pet_image",
        "t1_image_native",
        "mask_tissues",
        "flow_fields",
        "dartel_template",
        "reference_mask",
    ]

    # Unzipping
    for input_name, query in zip(["unzip_{_}" for _ in compressed_inputs], queries):
        wf.add(
            Nipype1Task(
                name=input_name,
                interface=Gunzip(),
                in_file=wf.lzin.T1w,  # TODO: change this
            )
        )

    # Coregister PET into T1 native space
    coreg_pet_t1 = Nipype1Task(
        name="coreg_pet_t1",
        interface=Coregister(),
    )
    coreg_pet_t1.inputs.in_files = wf.unzip_pet_image.lzout.out_file
    wf.add(coreg_pet_t1)

    # Spatially normalize PET into MNI
    dartel_mni_reg = Nipype1Task(
        name="dartel_mni_reg",
        interface=DARTELNorm2MNI(),
    )
    dartel_mni_reg.inputs.modulate = False
    dartel_mni_reg.inputs.fwhm = 0
    dartel_mni_reg.inputs.in_files = wf.unzip_flow_fields.lzout.out_file
    wf.add(dartel_mni_reg)

    # Reslice reference region mask into PET
    reslice = Nipype1Task(
        name="reslice",
        interface=Reslice(),
    )
    reslice.inputs.in_files = wf.unzip_reference_mask.lzout.out_file
    wf.add(reslice)

    # Normalize PET values according to reference region
    wf.add(
        normalize_to_reference(
            name="norm_to_ref",
            interface=normalize_to_reference,
            pet_image=wf.dartel_mni_reg.lzout.out_file,
            region_mask=wf.reslice.lzout.out_file,
        )
    )

    # Create binary mask from segmented tissues
    wf.add(
        create_binary_mask(
            name="binary_mask",
            interface=create_binary_mask,
            tissues=wf.unzip_mask_tissues.lzout.out_file,
            threshold=parameters["mask_threshold"],
        )
    )

    # Mask PET image
    wf.add(
        apply_binary_mask(
            name="apply_mask",
            interface=apply_binary_mask,
            image=wf.norm_to_ref.lzout.suvr_pet_path,
            binary_mask=wf.binary_mask.lzout.out_mask,
        )
    )

    # Smoothing
    if parameters["smooth"] is not None and len(parameters["smooth"]) > 0:
        smoothing_node = Nipype1Task(
            name="smoothing_node",
            interface=Smooth(),
        )
        smoothing_node.inputs.fwhm = [[x] * 3 for x in parameters["smooth"]]
        smoothing_node.inputs.out_prefix = [
            f"fwhm-{x}mm_" for x in parameters["smooth"]
        ]
        smoothing_node.inputs.in_file = wf.apply_mask.lzout.masked_image_path
        wf.add(smoothing_node)

    # Atlas Statistics
    wf.add(
        atlas_statistics(
            name="atlas_stats_node",
            interface=atlas_statistics,
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
        )

        # Creating Mask to use in PVC
        wf.add(
            create_pvc_mask(
                name="pvc_mask",
                interface=create_pvc_mask,
                tissues=wf.unzip_pvc_mask_tissues.lzout.out_file,
            )
        )

        # PET PVC
        petpvc = Nipype1Task(
            name="pvc",
            interface=PETPVC(),
        )
        petpvc.inputs.pvc = "RBV"
        petpvc.inputs.out_file = "pvc.nii"
        petpvc.inputs.in_file = wf.coreg_pet_t1.lzout.coregistered_source
        petpvc.inputs.mask_file = wf.pvc_mask.lzout.out_mask
        wf.add(petpvc)

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
        wf.add(dartel_mni_reg_pvc)

        # Reslice reference region mask into PET
        reslice_pvc = Nipype1Task(
            name="reslice_pvc",
            interfaces=Reslice(),
        )
        reslice_pvc.inputs.in_files = wf.unzip_reference_mask.lzout.out_file
        wf.add(reslice_pvc)

        # Normalize PET values according to reference region
        wf.add(
            normalize_to_reference(
                name="norm_to_ref_pvc",
                interface=normalize_to_reference,
                pet_image=wf.dartel_mni_reg_pvc.lzout.normalized_files,
                region_mask=wf.reslice_pvc.lzout.out_file,
            )
        )

        # Mask PET image
        wf.add(
            apply_binary_mask(
                name="apply_mask_pvc",
                interface=apply_binary_mask,
                image=wf.norm_to_ref_pvc.lzout.suvr_pet_path,
                binary_mask=wf.binary_mask.lzout.out_mask,
            )
        )

        # Smoothing
        if parameters["smooth"] is not None and len(parameters["smooth"]) > 0:
            smoothing_pvc = Nipype1Task(
                name="smoothing_pvc",
                interface=Smooth(),
            )
            smoothing_pvc.inputs.fwhm = [[x] * 3 for x in parameters["smooth"]]
            smoothing_pvc.inputs.out_prefix = [
                f"fwhm-{x}mm_" for x in parameters["smooth"]
            ]
            smoothing_pvc.inputs.in_file = wf.apply_mask_pvc.lzout.masked_image_path
            wf.add(smoothing_pvc)

        # Atlas Statistics
        wf.add(
            atlas_statistics(
                name="atlas_stats_pvc",
                interface=atlas_statistics,
                in_image=wf.norm_to_ref_pvc.lzout.suvr_pet_path,
                in_atlas_list=parameters["atlases"],
            )
        )

    wf.set_output(
        [
            ("pet_suvr_masked_smoothed", wf.smoothing_node.lzout.smoothed_files),
            ("pet_t1_native", wf.coreg_pet_t1.lzout.coregistered_source),
            ("pet_mni", wf.dartel_mni_reg.lzout.normalized_files),
            ("pet_suvr", wf.norm_to_ref.lzout.suvr_pet_path),
            ("binary_mask", wf.binary_mask.lzout.out_mask),
            ("pet_suvr_masked", wf.apply_mask.lzout.masked_image_path),
            ("atlas_statistics", wf.atlas_stats_node.lzout.atlas_statistics),
            ("pet_pvc_suvr_masked_smoothed", wf.smoothing_pvc.lzout.smoothed_files),
            ("pet_pvc", wf.petpvc.lzout.out_file),
            ("pet_pvc_mni", wf.dartel_mni_reg_pvc.lzout.normalized_files),
            ("pet_pvc_suvr", wf.norm_to_ref_pvc.lzout.suvr_pet_path),
            ("pet_pvc_suvr_masked", wf.apply_mask_pvc.lzout.masked_image_path),
            ("pvc_atlas_statistics", wf.atlas_stats_pvc.lzout.atlas_statistics),
        ]
    )

    return wf
