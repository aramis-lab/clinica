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


def run(wf: Workflow) -> None:

    """Execute a Pydra workflow

    Parameters
    ----------
        wf : Workflow
            The workflow to execute
    """

    with Submitter(plugin="cf") as submitter:
        submitter(wf)

    results = wf.result(return_inputs=True)
    # @TODO: decide where to store results
    print(results)
