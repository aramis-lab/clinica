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
    parameters = _check_pipeline_parameters(parameters)
    input_spec = pydra.specs.SpecInfo(
        name="Input",
        fields=[
            ("_graph_checksums", Any),
        ],
        bases=(pydra.specs.BaseSpec,),
    )
    wf = Workflow(
        name,
        input_spec=input_spec,
    )

    return wf
