from typing import Any

import pydra
from nipype.algorithms.misc import Gunzip
from pydra import Workflow
from pydra.tasks.nipype1.utils import Nipype1Task

from clinica.pydra.engine import clinica_io


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
    from clinica.utils.group import check_group_label

    if "group_label" not in parameters.keys():
        raise KeyError("Missing compulsory group_label key in pipeline parameter.")
    parameters.setdefault("dartel_tissues", [1, 2, 3])

    check_group_label(parameters["group_label"])

    return parameters


@clinica_io
def t1volume_register_dartel(
    name: str = "t1-volume-register-dartel", parameters: dict = {}
) -> Workflow:
    """Core workflow  for inter-subject registration using Dartel

    Parameters
    ----------
    name : str
        name of pipeline. Default="t1-volume-tissue-segmentation

    parameters : dict
        dictionary of pipeline parameters

    Returns
    -------
    Workflow
        pydra workflow for core functionalities

    """
    import clinica.pipelines.t1_volume_register_dartel.t1_volume_register_dartel_utils as utils
    from clinica.pydra.t1_volume.tasks import (
        task_create_iteration_parameters,
        task_prepare_dartel_input_images,
    )

    parameters = _check_pipeline_parameters(parameters)

    input_spec = pydra.specs.SpecInfo(
        name="Input",
        fields=[
            ("_graph_checksums", Any),
            (
                "dartel_input_tissue",
                dict,
                {"tissue_number": parameters["dartel_tissues"]},
                {"mandatory": True},
            ),
            (
                "dartel_iteration_templates",
                dict,
                {
                    "group_label": parameters["group_label"],
                    "i": range(1, 7),
                },
                {"mandatory": True},
            ),
        ],
        bases=(pydra.specs.BaseSpec,),
    )

    workflow = Workflow(
        name,
        input_spec=input_spec,
    )

    workflow.split("dartel_input_tissue")

    workflow.add(
        Nipype1Task(
            name="unzip_dartel_input_tissue",
            interface=Gunzip(),
            in_file=workflow.lzin.dartel_input_tissue,
        )
        .split("in_file")
        .combine("in_file")
    )

    workflow.add(
        Nipype1Task(
            name="unzip_templates",
            interface=Gunzip(),
            in_file=workflow.lzin.dartel_iteration_templates,
        )
        .split("in_file")
        .combine("in_file")
    )

    workflow.add(
        task_prepare_dartel_input_images(
            name="task_prepare_dartel_input_images",
            interface=task_prepare_dartel_input_images,
            nifti_input=workflow.unzip_dartel_input_tissue.lzout.out_file,
        )
    )

    workflow.add(
        task_create_iteration_parameters(
            name="task_create_iteration_parameters",
            template_input=workflow.unzip_templates.lzout.out_file,
        )
    )

    task_dartel_existing_template = Nipype1Task(
        name="task_dartel_existing_template",
        interface=utils.DARTELExistingTemplate(),
    )

    task_dartel_existing_template.inputs.image_files = (
        workflow.task_prepare_dartel_input_images.lzout.image_files
    )
    task_dartel_existing_template.inputs.iteration_parameters = (
        workflow.task_create_iteration_parameters.lzout.iteration_parameters
    )

    workflow.add(task_dartel_existing_template)

    workflow.set_output(
        [
            (
                "dartel_flow_fields",
                workflow.task_dartel_existing_template.lzout.dartel_flow_fields,
            )
        ]
    )

    return workflow
