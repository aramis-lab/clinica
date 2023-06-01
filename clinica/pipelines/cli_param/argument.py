"""Common CLI arguments used by Clinica pipelines."""
import click

from clinica.utils.pet import LIST_SUVR_REFERENCE_REGIONS, Tracer

acq_label = click.argument(
    "acq_label",
    type=click.Choice(Tracer),
)

bids_directory = click.argument(
    "bids_directory",
    type=click.Path(exists=True, file_okay=False, resolve_path=True),
)

caps_directory = click.argument(
    "caps_directory",
    type=click.Path(writable=True, file_okay=False, resolve_path=True),
)

contrast = click.argument("contrast")

group_label = click.argument("group_label")

orig_input_data = click.argument(
    "orig_input_data",
    type=click.Choice(["t1-freesurfer", "pet-surface", "custom-pipeline"]),
)

orig_input_data_volume = click.argument(
    "orig_input_data_volume",
    type=click.Choice(["t1-volume", "pet-volume", "custom-pipeline"]),
)

orig_input_data_ml = click.argument(
    "orig_input_data_ml",
    type=click.Choice(["t1-volume", "pet-surface"]),
)

pvc_psf_tsv = click.argument(
    "pvc_psf_tsv",
    type=click.Path(exists=True, dir_okay=False, resolve_path=True),
)

subject_visits_with_covariates_tsv = click.argument(
    "subject_visits_with_covariates_tsv",
    type=click.Path(exists=True, dir_okay=False, resolve_path=True),
)

suvr_reference_region = click.argument(
    "suvr_reference_region", type=click.Choice(LIST_SUVR_REFERENCE_REGIONS)
)
