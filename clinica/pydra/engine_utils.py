from pydra import Submitter, Workflow


def list_keys(query_bids: dict) -> list:
    """
    :query_bids: query dictionary to BIDSDataGrabber
    :return: IDS of queried elements
    """

    return list(query_bids.keys())


def list_out_fields(wf: Workflow) -> list:
    """
    :return: list of workflow "wf" fields in output_spec
    """
    return [x[0] for x in wf.output_spec.fields if not x[0].startswith("_")]


def list_in_fields(wf: Workflow) -> list:
    """
    :return: list of worfklow wf fields in input_spec
    """
    return [x[0] for x in wf.input_spec.fields if not x[0].startswith("_")]


def bids_query(keys: list) -> dict:
    """
    keys: a list of BIDS items to query
    returns: a dictionary formated for input to BIDSDataGrabber()
    """
    bids_keys_available = {
        "T1w": {"datatype": "anat", "suffix": "T1w", "extension": [".nii.gz"]}
    }

    return {key: bids_keys_available[key] for key in keys}


def run(wf):
    with Submitter(plugin="cf") as submitter:
        submitter(wf)

    results = wf.result(return_inputs=True)
    # @TODO: decide where to store results
    print(results)
