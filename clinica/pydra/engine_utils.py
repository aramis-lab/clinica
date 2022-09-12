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


def list_dict_in_fields(wf: Workflow) -> dict:
    """Extract the input default dictionary values from a Workflow.

    This is used as a way to pass data from the core Workflow
    to the input_workflow in order to properly query input folders.

    Better ways to do that ???

    Parameters
    ----------
    wf: Workflow

    Returns
    -------
    dict : Dictionary of wf default values from input_spec.
    """
    return {
        x[0]: x[2]
        for x in wf.input_spec.fields
        if not x[0].startswith("_") and x[1] == dict
    }


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

    return {key: bids_keys_available[key] for key in keys if key in bids_keys_available}


def caps_query(query: dict) -> dict:
    """Form dictionary of CAPS query based on a raw query.

    Parameters
    ----------
    query : Raw query dictionary
        The CAPS items to query as keys and the associated
        parameters as values.

    Returns
    -------
    dict
        Query dictionary compatible with CAPSDataGrabber()
    """
    from clinica.utils.input_files import (
        t1_volume_dartel_input_tissue,
        t1_volume_deformation_to_template,
        t1_volume_final_group_template,
        t1_volume_i_th_iteration_group_template,
        t1_volume_native_tpm,
        t1_volume_native_tpm_in_mni,
    )

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

    query_dict = {}
    for k, v in query.items():
        query_dict[k] = {}
        if k in caps_keys_available_file_reader:
            query_dict[k]["query"] = caps_keys_available_file_reader[k](**v)
            query_dict[k]["reader"] = "file"
        elif k in caps_keys_available_group_reader:
            query_dict[k]["query"] = caps_keys_available_group_reader[k](**v)
            query_dict[k]["reader"] = "group"
    return query_dict


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
