import typing as ty

import pydra
from nipype.interfaces.ants import ApplyTransforms, Registration, RegistrationSynQuick
from pydra.engine import Workflow
from pydra.tasks.nipype1.utils import Nipype1Task

from clinica.pydra.engine import clinica_io
from clinica.pydra.pet_linear.tasks import (
    concatenate_transforms_task,
    crop_nifti_task,
    suvr_normalization_task,
)
from clinica.pydra.tasks import download_mni_template, download_ref_template
from clinica.utils.pet import get_suvr_mask

IMAGE_DIMENSION = 3


@clinica_io
def build_core_workflow(name: str = "core", parameters: dict = {}) -> Workflow:
    """Build the core workflow for the PET Linear pipeline.

    Parameters
    ----------
    name : str, optional
        The name of the workflow. Default="core".

    parameters : dict, optional
        Optional dictionary of parameters to be used within the workflow.
        Default={}.

    Returns
    -------
    wf : Workflow
        The core workflow.
    """
    input_spec = pydra.specs.SpecInfo(
        name="Input",
        fields=[
            ("_graph_checksums", ty.Any),
            ("T1w", str, {"mandatory": True}),
            ("pet", dict, {"tracer": parameters["acq_label"]}, {"mandatory": True}),
            ("t1w_to_mni", dict, {}, {"mandatory": True}),
        ],
        bases=(pydra.specs.BaseSpec,),
    )

    ref_mask = get_suvr_mask(parameters["suvr_reference_region"])

    wf = Workflow(
        name,
        input_spec=input_spec,
    )

    wf.split(("pet", "T1w", "t1w_to_mni"))

    wf.add(download_mni_template(name="download_mni_template"))

    wf.add(download_ref_template(name="download_ref_template"))

    # RegistrationSynQuick by ANTS.
    ants_registration = Nipype1Task(
        name="ants_registration",
        interface=RegistrationSynQuick(),
    )
    ants_registration.inputs.dimension = IMAGE_DIMENSION
    ants_registration.inputs.transform_type = "r"
    ants_registration.inputs.fixed_image = wf.lzin.T1w
    ants_registration.inputs.moving_image = wf.lzin.pet
    wf.add(ants_registration)

    # Concatenate
    wf.add(
        concatenate_transforms_task(
            name="concatenate_transforms",
            interface=concatenate_transforms_task,
            pet_to_t1w_tranform=wf.ants_registration.lzout.out_matrix,
            t1w_to_mni_tranform=wf.lzin.t1w_to_mni,
        )
    )

    # ApplyTransforms by ANTS. PET to MRI
    apply_transform_pet_to_mni = Nipype1Task(
        name="apply_transform_pet_to_mni",
        interface=ApplyTransforms(),
    )
    apply_transform_pet_to_mni.inputs.dimension = IMAGE_DIMENSION
    apply_transform_pet_to_mni.inputs.reference_image = (
        wf.download_mni_template.lzout.mni_template_file
    )
    apply_transform_pet_to_mni.inputs.input_image = wf.lzin.pet
    apply_transform_pet_to_mni.inputs.transforms = (
        wf.concatenate_transforms.lzout.transforms_list
    )
    wf.add(apply_transform_pet_to_mni)

    # Normalize the image
    ants_registration_t1w_to_mni = Nipype1Task(
        name="ants_registration_t1w_to_mni",
        interface=Registration(),
    )
    ants_registration_t1w_to_mni.inputs.metric = ["MI"]
    ants_registration_t1w_to_mni.inputs.metric_weight = [1.0]
    ants_registration_t1w_to_mni.inputs.transforms = ["SyN"]
    ants_registration_t1w_to_mni.inputs.transform_parameters = [(0.1, 3, 0)]
    ants_registration_t1w_to_mni.inputs.dimension = IMAGE_DIMENSION
    ants_registration_t1w_to_mni.inputs.shrink_factors = [[8, 4, 2]]
    ants_registration_t1w_to_mni.inputs.smoothing_sigmas = [[3, 2, 1]]
    ants_registration_t1w_to_mni.inputs.sigma_units = ["vox"]
    ants_registration_t1w_to_mni.inputs.number_of_iterations = [[200, 50, 10]]
    ants_registration_t1w_to_mni.inputs.convergence_threshold = [1e-05]
    ants_registration_t1w_to_mni.inputs.convergence_window_size = [10]
    ants_registration_t1w_to_mni.inputs.radius_or_number_of_bins = [32]
    ants_registration_t1w_to_mni.inputs.winsorize_lower_quantile = 0.005
    ants_registration_t1w_to_mni.inputs.winsorize_upper_quantile = 0.995
    ants_registration_t1w_to_mni.inputs.collapse_output_transforms = True
    ants_registration_t1w_to_mni.inputs.use_histogram_matching = False
    ants_registration_t1w_to_mni.inputs.verbose = True
    ants_registration_t1w_to_mni.inputs.fixed_image = (
        wf.download_mni_template.lzout.mni_template_file
    )
    ants_registration_t1w_to_mni.inputs.moving_image = wf.lzin.T1w
    wf.add(ants_registration_t1w_to_mni)

    apply_transform_non_linear = Nipype1Task(
        name="apply_transform_non_linear",
        interface=ApplyTransforms(),
    )
    apply_transform_non_linear.inputs.dimension = IMAGE_DIMENSION
    apply_transform_non_linear.inputs.reference_image = (
        wf.download_mni_template.lzout.mni_template_file
    )
    apply_transform_non_linear.inputs.transforms = (
        wf.ants_registration_t1w_to_mni.lzout.reverse_forward_transforms
    )
    apply_transform_non_linear.inputs.input_image = (
        wf.apply_transform_pet_to_mni.lzout.output_image
    )
    wf.add(apply_transform_non_linear)

    wf.add(
        suvr_normalization_task(
            name="suvr_normalization",
            interface=suvr_normalization_task,
            input_img=wf.apply_transform_pet_to_mni.lzout.output_image,
            norm_img=wf.apply_transform_non_linear.lzout.output_image,
            ref_mask=ref_mask,
        )
    )

    output_connections = [
        ("affine_mat", wf.ants_registration.lzout.out_matrix),
        ("suvr_pet", wf.suvr_normalization.lzout.output_img),
    ]

    # Crop image
    if not parameters["uncropped_image"]:
        wf.add(
            crop_nifti_task(
                name="crop_nifti",
                interface=crop_nifti_task,
                input_img=wf.suvr_normalization.lzout.output_img,
                ref_img=wf.download_ref_template.lzout.ref_template_file,
            )
        )

        output_connections += [
            ("outfile_crop", wf.crop_nifti.lzout.output_img),
        ]

    # Optional argument
    if parameters["save_PETinT1w"]:
        apply_transform_pet_to_t1w = Nipype1Task(
            name="apply_transform_pet_to_t1w",
            interface=ApplyTransforms(),
        )
        apply_transform_pet_to_t1w.inputs.dimension = IMAGE_DIMENSION
        apply_transform_pet_to_t1w.inputs.input_image = wf.lzin.pet
        apply_transform_pet_to_t1w.inputs.reference_image = wf.lzin.T1w
        apply_transform_pet_to_t1w.inputs.transforms = (
            wf.ants_registration.lzout.out_matrix
        )

        wf.add(apply_transform_pet_to_t1w)

        output_connections += [
            ("PETinT1w", wf.apply_transform_pet_to_t1w.lzout.output_image),
        ]

    wf.set_output(output_connections)

    return wf
