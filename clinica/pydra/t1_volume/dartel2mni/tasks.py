from typing import Iterable, Union

from pydra.mark import annotate, task


@task
@annotate({"return": {"prepared_flowfields": list}})
def prepare_flowfields_task(flowfields: Iterable, number_of_tissues: int) -> list:
    """Pydra task for prepare_flow_fields."""
    return prepare_flowfields(flowfields, number_of_tissues)


def prepare_flowfields(flowfields: Iterable, number_of_tissues: int) -> list:
    """Reshape the flow fields according to the number of tissues.

    This function broadcasts the flow fields to the number of tissues
    selected by the user.

    Parameters
    ----------
    flowfields : list of str or str
        The path to the flowfields. It can be a string if we are dealing
        with a path to a single file, or a list if we are dealing with a
        list of paths.

    number_of_tissues : int
        The number of tissues selected.

    Returns
    -------
    list:
        Either a simple list with the flowfields path repeated for each tissue,
        of a list of lists of shape (n_flowfields, n_tissues).

    Raises
    ------
    TypeError
        If the provided flowfields is not a string, a list, or a tuple.
        If the provided number of tissues is not an integer.

    Examples
    --------
    >>> prepare_flowfields(["path_to_flowfield_1", "path_to_flowfield_2"], 3)
    [['path_to_flowfield_1', 'path_to_flowfield_1', 'path_to_flowfield_1'], ['path_to_flowfield_2', 'path_to_flowfield_2', 'path_to_flowfield_2']]
    >>> prepare_flowfields(["path_to_flowfield_1", "path_to_flowfield_2"], 2)
    [['path_to_flowfield_1', 'path_to_flowfield_1'], ['path_to_flowfield_2', 'path_to_flowfield_2']]
    >>> prepare_flowfields(["path_to_flowfield_1"], 2)
    [['path_to_flowfield_1', 'path_to_flowfield_1']]
    >>> prepare_flowfields("path_to_flowfield_1", 2)
    ['path_to_flowfield_1', 'path_to_flowfield_1']
    """
    if not isinstance(number_of_tissues, int):
        raise TypeError(
            f"Invalid type for number_of_tissues. Expected int, got {type(number_of_tissues)}."
        )
    if isinstance(flowfields, (list, tuple)):
        return [[f] * number_of_tissues for f in flowfields]
    if isinstance(flowfields, str):
        return [flowfields] * number_of_tissues
    raise TypeError(f"Invalid flowfields type: {type(flowfields)}.")
