from pathlib import Path
from typing import Iterable, List, Optional, Union

import click

bids_directory = click.argument(
    "bids_directory", type=click.Path(exists=True, resolve_path=True)
)

current_directory = click.argument(
    "current_directory", type=click.Path(exists=True, resolve_path=True)
)


@click.group("iotools")
def cli() -> None:
    """Tools to handle BIDS/CAPS datasets."""
    pass


@cli.command()
@click.argument("dataset", type=click.Path(resolve_path=True))
def describe(dataset: Union[str, Path]):
    """Describe a dataset in BIDS or CAPS format."""
    from .describe import describe as _describe

    _describe(dataset)


@cli.command()
@bids_directory
@click.argument("output_bids_directory", type=click.Path(resolve_path=True))
@click.option(
    "-m",
    "--modality",
    "modalities",
    multiple=True,
    metavar="modalities",
    help="Process selected modalities.",
    default=("T1w",),
)
@click.option("--center-all-files", is_flag=True, help="Force processing of all files.")
def center_nifti(
    bids_directory: str,
    output_bids_directory: str,
    modalities: Optional[Iterable[str]] = None,
    center_all_files: bool = False,
) -> None:
    """Center NIfTI files in a BIDS dataset."""
    import sys

    from clinica.utils.exceptions import ClinicaExistingDatasetError

    from .center_nifti import center_nifti as center_nifti_

    try:
        center_nifti_(
            bids_directory,
            output_bids_directory,
            modalities=modalities,
            center_all_files=center_all_files,
            overwrite_existing_files=False,  # todo : technically can remove this one
        )
    except ClinicaExistingDatasetError:
        click.echo(
            "Target BIDS directory is not empty. Existing files may be overwritten."
        )
        if click.confirm("Do you wish to continue?"):
            center_nifti_(
                bids_directory,
                output_bids_directory,
                modalities=modalities,
                center_all_files=center_all_files,
                overwrite_existing_files=True,
            )
        else:
            click.echo("Clinica will now exit...")
            sys.exit(0)


@cli.command()
@bids_directory
@click.argument("output_directory", type=click.Path(writable=True))
@click.option(
    "-op",
    "--output_prefix",
    default="missing_mods",
    show_default=True,
    help="Prefix for the name of the output files.",
)
def check_missing_modalities(
    bids_directory: str,
    output_directory: str,
    output_prefix: str = "missing_mods",
) -> None:
    """Check missing modalities in a BIDS dataset."""
    from clinica.iotools.utils.data_handling import compute_missing_mods
    from clinica.utils.inputs import check_bids_folder

    check_bids_folder(bids_directory)
    compute_missing_mods(bids_directory, output_directory, output_prefix)


@cli.command()
@bids_directory
@click.argument("caps_directory", type=click.Path(exists=True, resolve_path=True))
@click.argument("output_file", type=click.Path(resolve_path=True))
def check_missing_processing(
    bids_directory: str,
    caps_directory: str,
    output_file: str,
) -> None:
    """Check missing processing in a CAPS dataset."""
    from clinica.iotools.utils.data_handling import compute_missing_processing
    from clinica.utils.inputs import check_caps_folder

    check_caps_folder(caps_directory)
    compute_missing_processing(bids_directory, caps_directory, output_file)


@cli.command()
@current_directory
@click.argument("output_tsv", type=click.Path(resolve_path=True))
def create_subjects_visits(current_directory: str, output_tsv: str) -> None:
    """Export participants with their sessions."""
    from os import makedirs
    from os.path import basename, dirname

    from clinica.iotools.utils.data_handling import create_subs_sess_list
    from clinica.utils.inputs import determine_caps_or_bids
    from clinica.utils.stream import cprint

    is_bids = determine_caps_or_bids(current_directory)
    output_directory = dirname(output_tsv)
    makedirs(output_directory, exist_ok=True)
    create_subs_sess_list(
        current_directory, output_directory, basename(output_tsv), is_bids_dir=is_bids
    )
    cprint(f"The TSV file was saved to {output_tsv}.")


@cli.command()
@bids_directory
@click.argument("output_tsv", type=click.Path(resolve_path=True))
@click.option(
    "-caps",
    "--caps_directory",
    type=click.Path(exists=True),
    help="Path to a CAPS dataset.",
)
@click.option(
    "-p",
    "--pipeline",
    "pipelines",
    multiple=True,
    type=click.Choice(
        [
            "t1-freesurfer",
            "t1-volume",
            "pet-volume",
            "t1-freesurfer-longitudinal",
            "dwi-dti",
        ]
    ),
    help="Pipeline to merge to the output TSV file. All pipelines are merged by default.",
)
@click.option(
    "-vas",
    "--volume_atlas_selection",
    multiple=True,
    help="Atlas to merge for t1- and pet-volume. All atlases are merged by default.",
)
@click.option(
    "-fas",
    "--freesurfer_atlas_selection",
    multiple=True,
    help="Atlas to merge for t1-freesurfer. All atlases are merged by default",
)
@click.option(
    "-pvc",
    "--pvc_restriction",
    type=click.IntRange(0, 1),
    help=(
        "Restriction on the label [_pvc-rbv] for pet-volume. Default: Merge all atlases, "
        "0: Merge atlases without the label only, 1: Merge atlases with the label only."
    ),
)
@click.option(
    "-pts",
    "--pet_tracers_selection",
    multiple=True,
    help="PET tracer to merge. All PET tracers are merged by default.",
)
@click.option(
    "-group",
    "--group_selection",
    multiple=True,
    help="Group to merge. All groups are merged by default.",
)
@click.option(
    "-tsv",
    "--subjects_sessions_tsv",
    type=click.Path(exists=True, resolve_path=True),
    help="List of subjects and their sessions in TSV format.",
)
@click.option(
    "--ignore_scan_files",
    is_flag=True,
    help="Ignore scan files. This may accelerate the procedure.",
)
@click.option(
    "--ignore_session_scan_files",
    is_flag=True,
    help="Ignore session files. This may accelerate the procedure.",
)
def merge_tsv(
    bids_directory: str,
    output_tsv: str,
    caps_directory: Optional[str] = None,
    pipelines: Optional[List[str]] = None,
    volume_atlas_selection: Optional[List[str]] = None,
    freesurfer_atlas_selection: Optional[List[str]] = None,
    pvc_restriction: Optional[int] = None,
    pet_tracers_selection: Optional[List[str]] = None,
    group_selection: Optional[List[str]] = None,
    subjects_sessions_tsv: Optional[str] = None,
    ignore_scan_files: bool = False,
    ignore_session_scan_files: bool = False,
) -> None:
    """Merge clinical data into a single TSV file."""
    from .merge_tsv import merge_tsv as merge_tsv_

    merge_tsv_(
        bids_directory,
        output_tsv,
        caps_directory=caps_directory,
        pipelines=pipelines,
        volume_atlas_selection=volume_atlas_selection,
        freesurfer_atlas_selection=freesurfer_atlas_selection,
        pvc_restriction=pvc_restriction,
        pet_tracers_selection=pet_tracers_selection,
        group_selection=group_selection,
        subjects_sessions_tsv=subjects_sessions_tsv,
        ignore_scan_files=ignore_scan_files,
        ignore_session_scan_files=ignore_session_scan_files,
    )
