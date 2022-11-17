from os import PathLike

from pydra import Submitter, Workflow
from pydra.engine.specs import Result


def list_out_fields(wf: Workflow) -> list:
    """Extract the output fields from a Workflow

    Parameters
    ----------
    wf : Workflow

    Returns
    -------
    list
        list of wf workflow's fields from output_spec
    """
    return [x[0] for x in wf.output_spec.fields if not x[0].startswith("_")]


def list_workflow_inputs(wf: Workflow) -> dict:
    """Extract the input default dictionary values from a Workflow.

    This is used as a way to pass data from the core Workflow
    to the input_workflow in order to properly query input folders.

    Parameters
    ----------
    wf : Workflow
        Workflow to list inputs from.

    Returns
    -------
    dict : Dictionary of wf default values from input_spec.
    """
    inputs = {}
    for input_field in wf.input_spec.fields:
        if not input_field[0].startswith("_"):
            if input_field[1] == dict:
                inputs[input_field[0]] = input_field[2]
            elif input_field[1] == str:
                inputs[input_field[0]] = {}
            else:
                raise ValueError(
                    f"Could not parse input field {input_field}."
                    "Please review your input specifications."
                )
    return inputs


def run(wf: Workflow) -> Result:
    """Execute a Pydra workflow.

    If the execution of the workflow fails, the
    errors are printed on stdout.

    Parameters
    ----------
    wf : Workflow
        The workflow to be executed.

    Returns
    -------
    Result :
        The result of running the Workflow

    """
    import re

    try:
        with Submitter(plugin="cf") as submitter:
            submitter(wf)
    except Exception as e:
        path = re.search(r"/.*\.pklz", str(e))
        if path:
            print(read_error(path.group(0)))
        print(str(e))

    return wf.result(return_inputs=False)


def read_error(path: PathLike) -> str:
    """Read pklz file with error message

    Parameters
    ----------
    path : PathLike
        The path of the file containing the error message

    Returns
    -------
    str
        The error message
    """

    import cloudpickle as cp

    with open(path, "rb") as fp:
        err = cp.load(fp)
    return err["error message"][:-1]
