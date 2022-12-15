import pydra
from pydra.engine import Workflow

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
    from clinica.utils.atlas import T1_VOLUME_ATLASES
    from clinica.utils.group import check_group_label

    parameters.setdefault("group_label", None)
    check_group_label(parameters["group_label"])

    parameters.setdefault("atlases", T1_VOLUME_ATLASES)
    parameters.setdefault("modulate", True)

    return parameters


@clinica_io
def build_core_workflow(name: str = "core", parameters: dict = {}) -> Workflow:
    """Build the core workflow for the T1Volume-parcellation pipeline.

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
    from clinica.pydra.t1_volume.parcellation.tasks import atlas_statistics_task

    parameters = _check_pipeline_parameters(parameters)
    input_spec = pydra.specs.SpecInfo(
        name="Input",
        fields=[
            ("_graph_checksums", Any),
            (
                "t1_volume_template_tpm_in_mni",
                dict,
                {
                    "group_label": parameters["group_label"],
                    "tissue_number": 1,
                    "modulation": parameters["modulate"],
                },
                {"mandatory": True},
            ),
        ],
        bases=(pydra.specs.BaseSpec,),
    )
    wf = Workflow(
        name,
        input_spec=input_spec,
    )
    wf.add(
        atlas_statistics_task(
            name="atlas_statistics",
            interface=atlas_statistics_task,
            in_image=wf.lzin.file_list,
            atlas_list=parameters["atlases"],
        )
    )
    return wf
