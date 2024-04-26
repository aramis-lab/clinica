from __future__ import annotations

from os import PathLike
from pathlib import Path
from typing import Optional, Sequence

from pydra.engine import Workflow
from pydra.mark import annotate, task

__all__ = ["read_bids", "read_bids_dataset", "read_bids_files"]


@task
@annotate(
    {
        "return": {
            "dataset_description": Optional[dict],
            "participant_ids": list[str],
            "session_ids": Optional[list[str]],
        }
    }
)
def read_bids_dataset(dataset_path: PathLike):
    dataset_path = Path(dataset_path).resolve()
    description_file = dataset_path / "dataset_description.json"
    dataset_description = (
        json.load(description_file.read_text()) if description_file.exists() else None
    )

    try:
        _ = next(dataset_path.glob("*/ses-*"))
        multi_sessions = True
    except StopIteration:
        multi_sessions = False

    if multi_sessions:
        visits = dataset_path.glob("sub-*/ses-*")
        participant_ids, session_ids = list(
            map(
                list,
                zip(
                    *(
                        str(visit.relative_to(dataset_path)).split("/")
                        for visit in visits
                    )
                ),
            )
        )
    else:
        visits = dataset_path.glob("sub-*")
        participant_ids = sorted(
            str(visit.relative_to(dataset_path)) for visit in visits
        )
        session_ids = None

    return dataset_description, participant_ids, session_ids


@task
@annotate({"return": {"files": list[Path]}})
def read_bids_files(
    dataset_path: PathLike,
    participant_ids: Sequence[str] | None = None,
    session_ids: Sequence[str] | None = None,
    datatype: str | None = None,
    suffix: str | None = None,
    extension: str | None = None,
):
    dataset_path = Path(dataset_path).resolve()
    datatype = datatype or "*"
    suffix = suffix or "*"
    extension = extension or "*"
    files = []

    if all([participant_ids, session_ids]):
        for participant_id, session_id in zip([participant_ids, session_ids]):
            dir_pattern = f"{participant_id}/{session_id}/{datatype}"
            name_pattern = f"{participant_id}_{session_id}*_{suffix}.{extension}"
            file_pattern = f"{dir_pattern}/{name_pattern}"
            files += sorted(dataset_path.glob(file_pattern))
    elif participant_ids:
        for participant_id in participant_ids:
            dir_pattern = f"{participant_id}/**/{datatype}"
            name_pattern = f"{participant_id}*_{suffix}.{extension}"
            file_pattern = f"{dir_pattern}/{name_pattern}"
            files += sorted(dataset_path.glob(file_pattern))
    else:
        dir_pattern = f"**/{datatype}"
        name_pattern = f"*_{suffix}.{extension}"
        file_pattern = f"{dir_pattern}/{name_pattern}"
        files += sorted(dataset_path.glob(file_pattern))

    return files


def read_bids(output_queries: dict, **kwargs) -> Workflow:
    workflow = Workflow(name="read_bids", input_spec=["dataset_path"], **kwargs)

    workflow.add(
        read_bids_dataset(
            name="read_bids_dataset", dataset_path=workflow.lzin.dataset_path
        )
    )
    connections = {
        "dataset_description": workflow.read_bids_dataset.lzout.dataset_description,
        "participant_ids": workflow.read_bids_dataset.lzout.participant_ids,
        "session_ids": workflow.read_bids_dataset.lzout.session_ids,
    }

    for output_name, bids_query in output_queries.items():
        task_ = read_bids_files(
            name=f"read_{output_name}",
            dataset_path=workflow.lzin.dataset_path,
            participant_ids=workflow.read_bids_dataset.lzout.participant_ids,
            session_ids=workflow.read_bids_dataset.lzout.session_ids,
            **bids_query,
        )
        workflow.add(task_)
        connections.update({output_name: task_.lzout.files})

    workflow.set_output(connections=connections)

    return workflow
