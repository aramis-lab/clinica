import typing as ty
from os import PathLike
from pathlib import PurePath

import pydra
from nipype.interfaces.ants import N4BiasFieldCorrection, RegistrationSynQuick
from pydra import Workflow
from pydra.mark import annotate, task

from clinica.pydra.engine import clinica_io
from clinica.pydra.tasks import download_mni_template, download_ref_template

n4_bias_field_correction = N4BiasFieldCorrection(bspline_fitting_distance=300)
registration_syn_quick = RegistrationSynQuick(transform_type="a")


@task
@annotate({"return": {"cropped_image": PurePath}})
def crop_image(input_image: PathLike, template_image: PathLike) -> PurePath:
    from pathlib import Path

    from nilearn.image import resample_to_img

    cropped_image = Path(input_image).name.replace(".nii.gz", "_cropped.nii.gz")

    resample_to_img(
        source_img=str(input_image), target_img=str(template_image), force_resample=True
    ).to_filename(cropped_image)

    return cropped_image


@clinica_io
def build_core_workflow(name: str = "core", parameters={}) -> Workflow:
    """Core workflow for the T1 linear pipeline.

    Parameters
    ----------
    name : str
        The name of the workflow.

    Returns
    -------
    Workflow
        The core workflow.
    """
    from pydra.tasks.nipype1.utils import Nipype1Task

    input_spec = pydra.specs.SpecInfo(
        name="Input",
        fields=[("T1w", str, {"mandatory": True})],
        bases=(pydra.specs.BaseSpec,),
    )

    wf = Workflow(name, input_spec=input_spec)

    wf.add(download_mni_template(name="download_mni_template"))

    wf.add(download_ref_template(name="download_ref_template"))

    wf.add(
        Nipype1Task(
            name="n4_bias_field_correction",
            interface=n4_bias_field_correction,
            input_image=wf.lzin.T1w,
        ).split("input_image")
    )

    wf.add(
        Nipype1Task(
            name="registration_syn_quick",
            interface=registration_syn_quick,
            fixed_image=wf.download_mni_template.lzout.mni_template_file,
            moving_image=wf.n4_bias_field_correction.lzout.output_image,
        )
    )

    wf.add(
        crop_image(
            name="crop_image",
            interface=crop_image,
            input_image=wf.registration_syn_quick.lzout.warped_image,
            template_image=wf.download_ref_template.lzout.ref_template_file,
        )
    )

    wf.set_output(
        [
            ("t1w_corrected_file", wf.n4_bias_field_correction.lzout.output_image),
            ("t1w_registered_file", wf.registration_syn_quick.lzout.warped_image),
            ("t1w_cropped_file", wf.crop_image.lzout.cropped_image),
            ("xfm_file", wf.registration_syn_quick.lzout.out_matrix),
        ]
    )

    return wf
