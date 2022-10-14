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


def bids_query(raw_query: dict) -> dict:
    """Parse a raw BIDS query dictionary and return a properly
    formatted query dictionary that is compatible with the
    `BIDSDataGrabber` interface.

    Parameters
    ----------
    raw_query : dict
        The raw BIDS query as a dictionary. This may contain data
        that is not supported by the `BIDSDataGrabber`.

    Returns
    -------
    query : dict
        The formatted query dictionary compatible with `BIDSDataGrabber`.
    """
    query = {}
    bids_default_queries = {
        "T1w": {"datatype": "anat", "suffix": "T1w", "extension": [".nii.gz"]}
    }
    for k, q in raw_query.items():
        if k in bids_default_queries:
            query[k] = {**bids_default_queries[k], **q}
    return query


def caps_query(raw_query: dict) -> dict:
    """Parse a raw CAPS query dictionary and return a properly
    formatted query dictionary that is compatible with the
    `CAPSDataGrabber` interface.

    Parameters
    ----------
    raw_query : dict
        The raw CAPS query as a dictionary. This may contain data
        that is not supported by the `CAPSDataGrabber`.

    Returns
    -------
    query : dict
        The formatted query dictionary compatible with `CAPSDataGrabber`.
    """
    from clinica.utils.input_files import (
        t1_volume_dartel_input_tissue,
        t1_volume_deformation_to_template,
        t1_volume_final_group_template,
        t1_volume_i_th_iteration_group_template,
        t1_volume_native_tpm,
        t1_volume_native_tpm_in_mni,
    )

    query = {}
    caps_keys_available_file_reader = {
        "mask_tissues": t1_volume_native_tpm_in_mni,
        "flow_fields": t1_volume_deformation_to_template,
        "pvc_mask_tissues": t1_volume_native_tpm,
        "dartel_input_tissue": t1_volume_dartel_input_tissue,
    }
    caps_keys_available_group_reader = {
        "dartel_template": t1_volume_final_group_template,
        "dartel_iteration_templates": t1_volume_i_th_iteration_group_template,
    }
    for k, v in raw_query.items():
        if k in caps_keys_available_file_reader:
            query[k] = caps_keys_available_file_reader[k](**v)
            query[k]["reader"] = "file"
        elif k in caps_keys_available_group_reader:
            query[k] = caps_keys_available_group_reader[k](**v)
            query[k]["reader"] = "group"
    return query


def run(wf: Workflow) -> str:
    """Execute a Pydra workflow

    Parameters
    ----------
        wf : Workflow

    Returns
    -------
    str
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
