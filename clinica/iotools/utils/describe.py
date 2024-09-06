from pathlib import Path
from typing import Union

__all__ = ["describe"]


def describe(dataset: Union[str, Path]):
    """Describe a dataset in BIDS or CAPS format.

    Parameters
    ----------
    dataset : str or Path
        The path to the BIDS or CAPS dataset to describe.
    """
    import json

    from rich import box
    from rich.console import Console
    from rich.table import Table

    from clinica.utils.inputs import DatasetType

    dataset_description = Path(dataset) / "dataset_description.json"
    with open(dataset_description, "r") as fp:
        metadata = json.load(fp)

    is_caps = metadata["DatasetType"] == DatasetType.DERIVATIVE
    dataset_table = Table(title="Dataset information", box=box.ROUNDED)
    dataset_table.add_column("Name", justify="right", style="bright_cyan", no_wrap=True)
    dataset_table.add_column("Type", justify="right", style="bright_yellow")
    dataset_table.add_column("BIDS Specifications Version", style="bright_magenta")
    if is_caps:
        processing_metadata = metadata.pop("Processing")
        dataset_table.add_column(
            "CAPS Specifications Version", justify="right", style="bright_green"
        )
        dataset_table.add_row(
            metadata["Name"],
            metadata["DatasetType"],
            metadata["BIDSVersion"],
            metadata["CAPSVersion"],
        )
    else:
        dataset_table.add_row(
            metadata["Name"],
            metadata["DatasetType"],
            metadata["BIDSVersion"],
        )
    if is_caps:
        processing_table = Table(title="Processing information", box=box.ROUNDED)
        processing_table.add_column(
            "Name", justify="right", style="bright_cyan", no_wrap=True
        )
        processing_table.add_column("Date", style="bright_yellow")
        processing_table.add_column("Author", justify="right", style="bright_magenta")
        processing_table.add_column("Machine", justify="right", style="bright_green")
        processing_table.add_column("InputPath", justify="right", style="bright_red")
        processing_table.add_column(
            "ClinicaVersion", justify="right", style="bright_blue"
        )
        processing_table.add_column(
            "Dependencies", justify="right", style="bright_green"
        )

        for processing in processing_metadata:
            dependency_table = Table(title="Dependencies", box=box.ROUNDED)
            dependency_table.add_column("Name", justify="right")
            dependency_table.add_column("VersionConstraint", justify="right")
            dependency_table.add_column("InstalledVersion", justify="right")

            for dependency in processing["Dependencies"]:
                dependency_table.add_row(
                    dependency["Name"],
                    dependency["VersionConstraint"],
                    dependency["InstalledVersion"],
                )
            processing_table.add_row(
                processing["Name"],
                processing["Date"],
                processing["Author"],
                processing["Machine"],
                processing["InputPath"],
                processing["ClinicaVersion"],
                dependency_table,
            )
    console = Console()
    console.print(dataset_table)
    if is_caps:
        console.print(processing_table)
