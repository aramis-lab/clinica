from os import PathLike
from pathlib import Path
from typing import Optional

from pydra.mark import annotate, task


@task
@annotate(
    {
        "return": {
            "bids_file": Path,
            "participant_id": str,
            "session_id": Optional[str],
            "datatype": str,
            "suffix": str,
            "extension": str,
            "entities": dict,
        }
    }
)
def parse_bids_file(bids_file: PathLike):
    """Parse the components of a BIDS file.

    Examples
    --------
    >>> task = parse_bids_file(bids_file=Path("dataset/sub-P01/ses-M00/anat/sub-P01_ses-M00_desc-test_T1w.nii.gz"))
    >>> result = task()
    >>> result.output.participant_id
    'sub-P01'
    >>> result.output.session_id
    'ses-M00'
    >>> result.output.datatype
    'anat'
    >>> result.output.suffix
    'T1w'
    >>> result.output.extension
    'nii.gz'
    >>> result.output.entities
    {'desc': 'test'}
    """
    bids_file = Path(bids_file)

    datatype = bids_file.parent.name
    file_name = bids_file.name
    rest, extension = file_name.split(sep=".", maxsplit=1)
    rest, suffix = rest.rsplit(sep="_", maxsplit=1)
    entities = {
        prefix: value for prefix, value in (part.split("-") for part in rest.split("_"))
    }

    participant_id = f"sub-{entities.pop('sub')}"
    try:
        session_id = f"ses-{entities.pop('ses')}"
    except KeyError:
        session_id = None

    return bids_file, participant_id, session_id, datatype, suffix, extension, entities
