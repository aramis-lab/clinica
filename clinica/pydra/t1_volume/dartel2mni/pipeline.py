import pydra
from nipype.algorithms.misc import Gunzip
from nipype.interfaces.spm import DARTELNorm2MNI
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

    if "group_label" not in parameters:
        raise KeyError("Missing compulsory group_label key in pipeline parameter.")
    check_group_label(parameters["group_label"])
    parameters.setdefault("tissues", [1, 2, 3])
    parameters.setdefault("voxel_size", None)
    parameters.setdefault("modulate", True)
    parameters.setdefault("smooth", 8.0)

    return parameters


@clinica_io
def t1volume_dartel2mni(
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
    from clinica.utils.spm import spm_standalone_is_available, use_spm_standalone

    if spm_standalone_is_available():
        use_spm_standalone()

    parameters = _check_pipeline_parameters(parameters)

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
    wf.split("in_file")
    compressed_inputs = [
        "tissues",
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

    # DARTEL2MNI Registration
    dartel2mni_node = Nipype1Task(
        name="dartel_2mni",
        interface=DARTELNorm2MNI(),
    )
    if parameters["voxel_size"] is not None:
        dartel2mni_node.inputs.voxel_size = tuple(parameters["voxel_size"])
    dartel2mni_node.inputs.modulate = parameters["modulate"]
    dartel2mni_node.inputs.fwhm = 0
    dartel2mni_node.inputs.apply_to_files = wf.unzip_native_segmentations.lzout.out_file

    # Smoothing
    if parameters["smooth"] is not None:
        from clinica.pydra.pet_volume.pipeline import (
            _sanitize_fwhm,
            build_smoothing_workflow,
        )

        smoothing_fwhm = _sanitize_fwhm(parameters["smooth"])
        smoothing_wf = build_smoothing_workflow(
            "smoothing_workflow",
            smoothing_fwhm,
        )
        smoothing_wf.inputs.input_file = wf.dartel2mni_node.lzout.normalized_files
        wf.add(smoothing_wf)
