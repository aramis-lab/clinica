from typing import Any

import pydra
from nipype.algorithms.misc import Gunzip
from nipype.interfaces.spm import DARTELNorm2MNI
from pydra import Workflow
from pydra.tasks.nipype1.utils import Nipype1Task

from clinica.pydra.engine import clinica_io


@clinica_io
def build_core_workflow(
    name: str = "t1-volume-dartel2mni", parameters: dict = {}
) -> Workflow:
    """Core workflow for T1VolumeDartel2MNI - Dartel template to MNI (Pydra engine).

    Parameters
    ----------
    name : str, optional
        Name of the workflow. Default="t1-volume-dartel2mni".

    parameters : dict, optional
        Dictionary of parameters. Default={}.
    """
    from clinica.pydra.shared_workflows.smoothing import build_smoothing_workflow
    from clinica.pydra.t1_volume.dartel2mni.tasks import prepare_flowfields_task
    from clinica.pydra.utils import sanitize_fwhm
    from clinica.utils.spm import spm_standalone_is_available, use_spm_standalone

    if spm_standalone_is_available():
        use_spm_standalone()

    input_spec = pydra.specs.SpecInfo(
        name="Input",
        fields=[
            ("_graph_checksums", Any),
            (
                "tissues",
                dict,
                {"tissue_number": parameters["tissues"]},
                {"mandatory": False},
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
        ],
        bases=(pydra.specs.BaseSpec,),
    )

    wf = Workflow(name, input_spec=input_spec)

    wf.split(("tissues", "flow_fields"))

    compressed_inputs = [
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
            name="unzip_tissues",
            interface=Gunzip(),
            in_file=wf.lzin.tissues,
        )
        .split("in_file")
        .combine("in_file")
    )

    # DARTEL2MNI Registration
    wf.add(
        prepare_flowfields_task(
            name="prepare_flowfields",
            interface=prepare_flowfields_task,
            flowfields=wf.unzip_flow_fields.lzout.out_file,
            number_of_tissues=len(parameters["tissues"]),
        )
    )

    dartel2mni_node = Nipype1Task(
        name="dartel2mni",
        interface=DARTELNorm2MNI(),
    )
    if parameters["voxel_size"] is not None:
        dartel2mni_node.inputs.voxel_size = tuple(parameters["voxel_size"])
    dartel2mni_node.inputs.modulate = parameters["modulate"]
    dartel2mni_node.inputs.fwhm = 0
    dartel2mni_node.inputs.apply_to_files = wf.unzip_tissues.lzout.out_file
    dartel2mni_node.inputs.template_file = wf.unzip_dartel_template.lzout.out_file
    dartel2mni_node.inputs.flowfield_files = (
        wf.prepare_flowfields.lzout.prepared_flowfields
    )
    wf.add(dartel2mni_node)

    output_connections = [("normalized_files", wf.dartel2mni.lzout.normalized_files)]

    # Smoothing
    if parameters["smooth"] is not None:
        smoothing_wf = build_smoothing_workflow(
            "smoothing_workflow",
            sanitize_fwhm(parameters["smooth"]),
        )
        smoothing_wf.inputs.input_file = wf.dartel2mni.lzout.normalized_files
        wf.add(smoothing_wf)
        output_connections += [
            ("smoothed_normalized_files", wf.smoothing_workflow.lzout.smoothed_files),
        ]

    wf.set_output(output_connections)

    return wf
