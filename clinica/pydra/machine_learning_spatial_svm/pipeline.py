import pydra
from pydra import Workflow

from clinica.pydra.engine import clinica_io


@clinica_io
def build_core_workflow(name: str = "core", parameters={}) -> Workflow:
    """Build the core workflow for the machine learning spatial svm pipeline.

    Parameters
    ----------
    name : str, optional
        The name of the workflow. Default="core".

    parameters : dict, optional
        Dictionary of parameters to be used
        within the workflow.
        Default={}.

    Returns
    -------
    wf : Workflow
        The core
    """
    from typing import Any

    input_spec = pydra.specs.SpecInfo(
        name="Input",
        fields=[
            ("_graph_checksums", Any),
            (query_name, dict, query, {"mandatory": True}),
        ],
        basqes=(pydra.specs.BaseSpec,),
    )
    wf = Workflow(name, input_spec=input_spec)
    return wf
