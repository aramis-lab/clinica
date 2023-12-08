from os import PathLike
from pathlib import PurePath
from typing import Optional

from pydra.mark import annotate, task


@task
@annotate(
    {
        "return": {
            "bids_file": PurePath,
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
    bids_file = PurePath(bids_file)

    datatype = bids_file.parent.name
    file_name = bids_file.name
    rest, extension = file_name.split(sep=".", maxsplit=1)
    rest, suffix = rest.rsplit(sep="_", maxsplit=1)
    entities = {
        prefix: value for prefix, value in (part.split("-") for part in rest.split("_"))
    }

    participant_id = entities.pop("sub")
    try:
        session_id = entities.pop("ses")
    except KeyError:
        session_id = None

    return bids_file, participant_id, session_id, datatype, suffix, extension, entities
