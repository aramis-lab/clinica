from pathlib import Path
from typing import Optional, Union

from rich import box
from rich.console import Console
from rich.table import Table

from clinica.converters.bids_dataset_description import (
    BIDSDatasetDescription,  # TODO: move that
)
from clinica.utils.caps import CAPSDatasetDescription
from clinica.utils.inputs import DatasetType

__all__ = ["describe"]


def describe(dataset: Union[str, Path]) -> None:
    """Describe a dataset in BIDS or CAPS format.

    Parameters
    ----------
    dataset : str or Path
        The path to the BIDS or CAPS dataset to describe.
    """
    description = (
        CAPSDatasetDescription if _is_caps(dataset) else BIDSDatasetDescription
    ).from_file(Path(dataset) / "dataset_description.json")
    console = Console()
    console.print(_build_dataset_table(description))
    if (processing_table := _build_processing_table(description)) is not None:
        console.print(processing_table)


def _is_caps(dataset: Union[str, Path]) -> bool:
    import json

    from clinica.utils.inputs import DatasetType

    dataset_description = Path(dataset) / "dataset_description.json"
    with open(dataset_description, "r") as fp:
        metadata = json.load(fp)

    return metadata["DatasetType"] == DatasetType.DERIVATIVE


def _build_dataset_table(
    description: Union[BIDSDatasetDescription, CAPSDatasetDescription],
) -> Table:
    dataset_table = Table(title="Dataset information", box=box.ROUNDED)
    for column, color in zip(
        ["Name", "Type", "BIDS Specifications Version"],
        ["bright_cyan", "bright_yellow", "bright_magenta"],
    ):
        dataset_table.add_column(column, justify="right", style=color)
    row = (description.name, description.dataset_type, str(description.bids_version))
    if description.dataset_type == DatasetType.DERIVATIVE:
        dataset_table.add_column(
            "CAPS Specifications Version", justify="right", style="bright_green"
        )
        row = (
            description.name,
            description.dataset_type,
            str(description.bids_version),
            str(description.caps_version),
        )
    dataset_table.add_row(*row)

    return dataset_table


def _build_processing_table(
    description: Union[BIDSDatasetDescription, CAPSDatasetDescription],
) -> Optional[Table]:
    if description.dataset_type == DatasetType.RAW:
        return None
    processing_table = Table(title="Processing information", box=box.ROUNDED)
    for column, color in zip(
        [
            "Name",
            "Date",
            "Author",
            "Machine",
            "InputPath",
            "ClinicaVersion",
            "Dependencies",
        ],
        [
            "bright_cyan",
            "bright_yellow",
            "bright_magenta",
            "bright_green",
            "bright_red",
            "bright_blue",
            "bright_green",
        ],
    ):
        processing_table.add_column(column, justify="right", style=color)
    for processing in description.processing:
        dependency_table = Table(title="Dependencies", box=box.ROUNDED)
        for column in ("Name", "VersionConstraint", "InstalledVersion"):
            dependency_table.add_column(column, justify="right")
        for dependency in processing.dependencies:
            dependency_table.add_row(
                dependency.name,
                str(dependency.version_constraint),
                str(dependency.installed_version),
            )
        processing_table.add_row(
            processing.name,
            str(processing.date),
            processing.author,
            processing.machine,
            processing.input_path,
            str(processing.clinica_version),
            dependency_table,
        )
    return processing_table
