from typing import Iterable

from pydra.mark import annotate, task


@task
@annotate({"return": {"prepared_flowfields": list}})
def prepare_flowfields_task(flowfields: Iterable, tissues: Iterable) -> list:
    """Pydra task for prepare_flow_fields."""
    return prepare_flowfields(flowfields, tissues)


def prepare_flowfields(flowfields: Iterable, tissues: Iterable) -> list:
    """Reshape the flow fields according to the number of tissues.

    This function broadcasts the flow fields to the number of tissues
    selected by the user.

    Parameters
    ----------
    flowfields : list of str or str
        The path to the flowfields. It can be a string if we are dealing
        with a path to a single file, or a list if we are dealing with a
        list of paths.

    tissues : Iterable
        The tissues selected. It should be a list of integers encoding the
        different tissues.

    Returns
    -------
    list:
        Either a simple list with the flowfields path repeated for each tissue,
        of a list of lists of shape (n_flowfields, n_tissues).

    Raises
    ------
    ValueError
        If the provided flowfields is not a string, a list, or a tuple.
        If the provided tissues is not an iterable.

    Examples
    --------
    >>> prepare_flowfields(["path_to_flowfield_1", "path_to_flowfield_2"], (1, 2, 3))
    [['path_to_flowfield_1', 'path_to_flowfield_1', 'path_to_flowfield_1'], ['path_to_flowfield_2', 'path_to_flowfield_2', 'path_to_flowfield_2']]
    >>> prepare_flowfields(["path_to_flowfield_1", "path_to_flowfield_2"], (1, 3))
    [['path_to_flowfield_1', 'path_to_flowfield_1'], ['path_to_flowfield_2', 'path_to_flowfield_2']]
    >>> prepare_flowfields(["path_to_flowfield_1"], (1, 3))
    [['path_to_flowfield_1', 'path_to_flowfield_1']]
    >>> prepare_flowfields("path_to_flowfield_1", (1, 3))
    ['path_to_flowfield_1', 'path_to_flowfield_1']
    """
    if not isinstance(tissues, Iterable):
        raise ValueError(f"Invalid tissues type: {type(tissues)}.")
    if isinstance(flowfields, (list, tuple)):
        return [[f] * len(tissues) for f in flowfields]
    if isinstance(flowfields, str):
        return [flowfields] * len(tissues)
    raise ValueError(f"Invalid flowfields type: {type(flowfields)}.")
