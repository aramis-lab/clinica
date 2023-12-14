import typing as ty
from os import PathLike
from pathlib import PurePath

import pydra
from nipype.interfaces.ants import N4BiasFieldCorrection, RegistrationSynQuick
from pydra import Workflow
from pydra.mark import annotate, task
from pydra.tasks.bids import read_bids_dataset

from clinica.pydra.engine import clinica_io
from clinica.pydra.tasks import (
    download_mni_template_2009c,
    download_ref_template,
    parse_bids_file,
    write_bids_file,
    write_caps_file,
)

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


# @clinica_io
def build_core_workflow(name: str = "core", **kwargs) -> Workflow:
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
        fields=[
            ("input_dir", str, {"mandatory": True}),
            ("output_dir", str, {"mandatory": True}),
        ],
        bases=(pydra.specs.BaseSpec,),
    )

    wf = Workflow(name=name, input_spec=input_spec, **kwargs)

    wf.add(download_mni_template_2009c(name="download_mni_template"))

    wf.add(download_ref_template(name="download_ref_template"))

    wf.add(
        read_bids_dataset(
            name="read_t1w",
            dataset_path=wf.lzin.input_dir,
            output_queries={"t1w_file": {"suffix": "T1w", "extension": ["nii", "nii.gz"]}},
        )
    )

    wf.add(parse_bids_file(name="parse_t1w", bids_file=wf.read_t1w.lzout.t1w_file).split("bids_file"))

    wf.add(
        Nipype1Task(
            name="n4_bias_field_correction",
            interface=n4_bias_field_correction,
            input_image=wf.parse_t1w.lzout.bids_file,
        )
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

    wf.add(
        write_bids_file(
            name="write_bias_corrected_t1w_file",
            input_file=wf.n4_bias_field_correction.lzout.output_image,
            dataset_path=wf.lzin.output_dir,
            participant_id=wf.parse_t1w.lzout.participant_id,
            session_id=wf.parse_t1w.lzout.session_id,
            datatype="anat",
            suffix="T1w",
            entities={"desc": "bfc"},
        )
    )

    wf.add(
        write_bids_file(
            name="write_registered_t1w_file",
            input_file=wf.registration_syn_quick.lzout.warped_image,
            dataset_path=wf.lzin.output_dir,
            participant_id=wf.parse_t1w.lzout.participant_id,
            session_id=wf.parse_t1w.lzout.session_id,
            datatype="anat",
            suffix="T1w",
            entities={"space": "MNI152NLin2009cSym", "res": "1x1x1", "desc": "nocrop"},
        )
    )

    wf.add(
        write_bids_file(
            name="write_cropped_t1w_file",
            input_file=wf.crop_image.lzout.cropped_image,
            dataset_path=wf.lzin.output_dir,
            participant_id=wf.parse_t1w.lzout.participant_id,
            session_id=wf.parse_t1w.lzout.session_id,
            datatype="anat",
            suffix="T1w",
            entities={"space": "MNI152NLin2009cSym", "res": "1x1x1", "desc": "cropped"},
        )
    )

    wf.add(
        write_caps_file(
            name="write_xfm_file",
            input_file=wf.registration_syn_quick.lzout.out_matrix,
            dataset_path=wf.lzin.output_dir,
            participant_id=wf.parse_t1w.lzout.participant_id,
            session_id=wf.parse_t1w.lzout.session_id,
            datatype="t1linear",
            suffix="xfm_space-{space}_affine",
            entities={"space": "MNI152NLin2009cSym"},
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
