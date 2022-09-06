from os import PathLike

from pydra import Submitter, Workflow


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


def list_in_fields(wf: Workflow) -> list:
    """Extract the input fields from a Workflow

    Parameters
    ----------
    wf : Workflow

    Returns
    -------
    list
        list of wf workflow's fields from input_spec
    """
    return [x[0] for x in wf.input_spec.fields if not x[0].startswith("_")]


def bids_query(keys: list) -> dict:
    """Form dictionary of BIDS query based on modality or key

    Parameters
    ----------
    keys : list
            The BIDS items to query
    Returns
    -------
    dict
            Query dictionary compatible with BIDSDataGrabber()
    """
    bids_keys_available = {
        "T1w": {"datatype": "anat", "suffix": "T1w", "extension": [".nii.gz"]}
    }

    return {key: bids_keys_available[key] for key in keys}


def run(wf: Workflow) -> str:
    """Execute a Pydra workflow

    Parameters
    ----------
        wf : Workflow

    Returns
    -------
    Union[str, pydra.engine.specs.Result]
        The result of running the Workflow

    """
    import re

    try:
        with Submitter(plugin="cf") as submitter:
            submitter(wf)
    except Exception as e:
        path = re.search("\/.*\.pklz", str(e))
        if path:
            return read_error(path.group(0))
        return str(e)

    results = wf.result(return_inputs=False)
    return str(results)


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
