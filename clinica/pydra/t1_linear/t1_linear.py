from os import PathLike
from pathlib import PurePath
from typing import Any

import pydra
from pydra import Workflow
from pydra.mark import annotate, task

from clinica.pydra.engine import clinica_io
from clinica.pydra.tasks import download_mni_template_2009c, download_ref_template


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
    from pydra.tasks import ants

    input_spec = pydra.specs.SpecInfo(
        name="Input",
        fields=[
            ("_graph_checksums", Any),
            ("T1w", str, {"mandatory": True}),
        ],
        bases=(pydra.specs.BaseSpec,),
    )

    wf = Workflow(name=name, input_spec=input_spec)

    wf.add(download_mni_template_2009c(name="download_mni_template"))

    wf.add(download_ref_template(name="download_ref_template"))

    wf.add(
        ants.N4BiasFieldCorrection(
            name="n4_bias_field_correction",
            dimensionality=3,
            input_image=wf.lzin.T1w,
            spline_distance=600,
        ).split("input_image")
    )

    wf.add(
        ants.registration_syn_quick(
            name="registration_syn_quick",
            dimensionality=3,
            transform_type="a",
            fixed_image=wf.download_mni_template.lzout.mni_template_file,
            moving_image=wf.n4_bias_field_correction.lzout.output_image,
        )
    )

    wf.add(
        ants.ApplyTransforms(
            name="apply_transforms",
            dimensionality=3,
            fixed_image=wf.download_ref_template.lzout.ref_template_file,
            moving_image=wf.registration_syn_quick.lzout.warped_moving_image,
            output_datatype="short",
        )
    )

    wf.set_output(
        [
            ("t1w_file", wf.apply_transforms.lzout.output_image),
            ("xfm_file", wf.registration_syn_quick.lzout.affine_transform),
        ]
    )

    return wf
