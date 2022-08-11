from typing import Optional

from pydra.engine.task import FunctionTask

DEFAULT_QUERY = {
    "bold": {
        "suffix": "bold",
        "extension": [".nii", ".nii.gz"],
    },
    "T1w": {
        "suffix": "T1w",
        "extension": [".nii", ".nii.gz"],
    },
}


def read_bids(output_query: Optional[dict] = None, **kwargs) -> FunctionTask:
    """Generate a BIDS reading task.

    Parameters
    ----------
    output_query : dict
        Mapping between output and BIDS query
    **kwargs
        Additional arguments passed to the task constructor

    Returns
    -------
    pydra.engine.task.FunctionTask
        BIDS reading task

    Examples
    --------
    >>> task = read_bids()
    >>> task.output_names
    ['bold', 'T1w']
    """
    from os import PathLike

    from pydra.engine.specs import BaseSpec, SpecInfo

    output_query = output_query or DEFAULT_QUERY

    def inner(base_dir):
        from ancpbids import BIDSLayout

        layout = BIDSLayout(ds_dir=base_dir)

        results = [
            layout.get(return_type="files", **query)
            for key, query in list(output_query.items())
        ]

        return tuple(results) if len(results) > 1 else results[0]

    input_spec = SpecInfo(
        name="Input",
        fields=[
            (
                "base_dir",
                PathLike,
                {"help_string": "Root directory of the BIDS dataset"},
            ),
        ],
        bases=(BaseSpec,),
    )

    output_spec = SpecInfo(
        name="Output",
        fields=[(key, str) for key in list(output_query.keys())],
        bases=(BaseSpec,),
    )

    return FunctionTask(
        func=inner,
        input_spec=input_spec,
        output_spec=output_spec,
        **kwargs,
    )
