from typing import List, Optional

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
@bids_directory
@click.argument("output_bids_directory", type=click.Path(resolve_path=True))
@click.option(
    "-m",
    "--modality",
    "modalities",
    multiple=True,
    metavar="modalities",
    help="Process selected modalities.",
)
@click.option("--center-all-files", is_flag=True, help="Force processing of all files.")
def center_nifti(
    bids_directory: str,
    output_bids_directory: str,
    modalities: Optional[List[str]] = None,
    center_all_files: bool = False,
) -> None:
    """Center NIfTI files in a BIDS dataset.

    This tool is mainly useful as a preprocessing step of SPM. In some cases, SPM is not able to segment T1 volumes
    because their respective center is not aligned with the origin of the world coordinate system. By default, only
    images detected as problematic are converted from INPUT_BIDS_DIRECTORY to OUTPUT_BIDS_DIRECTORY, whilst the others
    are copied verbatim.
    """
    import sys
    import time
    from os import listdir, makedirs
    from os.path import isfile, join

    from clinica.iotools.utils.data_handling import (
        center_all_nifti,
        write_list_of_files,
    )
    from clinica.utils.stream import cprint

    # check that output_folder does not exist, or is an empty folder
    try:
        makedirs(output_bids_directory)
    except FileExistsError:
        file_list = [
            file for file in listdir(output_bids_directory) if not file.startswith(".")
        ]
        if file_list:
            click.echo(
                "Target BIDS directory is not empty. Existing files may be overwritten."
            )
            if not click.confirm("Do you wish to continue?"):
                click.echo("Clinica will now exit...")
                sys.exit(0)

    cprint("Clinica is now centering the requested images.")

    centered_files = center_all_nifti(
        bids_directory,
        output_bids_directory,
        modalities,
        center_all_files,
    )

    # Write list of created files
    timestamp = time.strftime("%Y%m%d-%H%M%S", time.localtime(time.time()))
    log_file = join(output_bids_directory, "centered_nifti_list_" + timestamp + ".txt")
    # If an error happen while creating the file, the function returns Nan
    if not write_list_of_files(centered_files, log_file):
        cprint("Could not create log file.")

    # Final message
    cprint(
        f"{str(len(centered_files))} NIfTI files/images of BIDS folder:\n"
        f"\t{bids_directory}\n"
        f"for the modalities {modalities} have been centered in output folder:\n"
        f"\t{output_bids_directory}"
    )
    if isfile(log_file):
        cprint(f"The list of centered NIfTI files is available here: {log_file}.")
    cprint(
        "Please note that the rest of the input BIDS folder has also been copied to the output folder."
    )


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
    help="Pipeline to merge to the ouput TSV file. All pipelines are merged by default.",
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

    from clinica.iotools.utils.data_handling import create_merge_file
    from clinica.utils.inputs import check_bids_folder

    check_bids_folder(bids_directory)

    create_merge_file(
        bids_directory,
        output_tsv,
        caps_dir=caps_directory,
        pipelines=pipelines,
        ignore_scan_files=ignore_scan_files,
        ignore_sessions_files=ignore_session_scan_files,
        volume_atlas_selection=volume_atlas_selection,
        freesurfer_atlas_selection=freesurfer_atlas_selection,
        pvc_restriction=pvc_restriction,
        tsv_file=subjects_sessions_tsv,
        group_selection=group_selection,
        tracers_selection=pet_tracers_selection,
    )
