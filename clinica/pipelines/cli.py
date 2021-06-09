import click

from .deeplearning_prepare_data import deeplearning_prepare_data_cli
from .dwi_connectome import dwi_connectome_cli
from .dwi_dti import dwi_dti_cli
from .dwi_preprocessing_using_phasediff_fieldmap import (
    dwi_preprocessing_using_phasediff_fieldmap_cli,
)
from .dwi_preprocessing_using_t1 import dwi_preprocessing_using_t1_cli
from .machine_learning_spatial_svm import spatial_svm_cli
from .pet_linear import pet_linear_cli
from .pet_surface import pet_surface_cli
from .pet_volume import pet_volume_cli
from .statistics_surface import statistics_surface_cli
from .statistics_volume import statistics_volume_cli
from .statistics_volume_correction import statistics_volume_correction_cli
from .t1_freesurfer import t1_freesurfer_cli
from .t1_freesurfer_longitudinal import (
    t1_freesurfer_longitudinal_cli,
    t1_freesurfer_longitudinal_correction_cli,
    t1_freesurfer_template_cli,
)
from .t1_linear import t1_linear_cli
from .t1_volume import t1_volume_cli
from .t1_volume_create_dartel import t1_volume_create_dartel_cli
from .t1_volume_existing_template import t1_volume_existing_template_cli
from .t1_volume_parcellation import t1_volume_parcellation_cli
from .t1_volume_register_dartel import t1_volume_register_dartel_cli
from .t1_volume_tissue_segmentation import t1_volume_tissue_segmentation_cli


@click.group(name="run")
def cli() -> None:
    pass


cli.add_command(deeplearning_prepare_data_cli.cli)
cli.add_command(dwi_connectome_cli.cli)
cli.add_command(dwi_dti_cli.cli)
cli.add_command(dwi_preprocessing_using_phasediff_fieldmap_cli.cli)
cli.add_command(dwi_preprocessing_using_t1_cli.cli)
cli.add_command(spatial_svm_cli.cli)
cli.add_command(pet_linear_cli.cli)
cli.add_command(pet_surface_cli.cli)
cli.add_command(pet_volume_cli.cli)
cli.add_command(statistics_surface_cli.cli)
cli.add_command(statistics_volume_cli.cli)
cli.add_command(statistics_volume_correction_cli.cli)
cli.add_command(t1_freesurfer_cli.cli)
cli.add_command(t1_freesurfer_longitudinal_cli.cli)
cli.add_command(t1_freesurfer_longitudinal_correction_cli.cli)
cli.add_command(t1_freesurfer_template_cli.cli)
cli.add_command(t1_linear_cli.cli)
cli.add_command(t1_volume_cli.cli)
cli.add_command(t1_volume_create_dartel_cli.cli)
cli.add_command(t1_volume_existing_template_cli.cli)
cli.add_command(t1_volume_parcellation_cli.cli)
cli.add_command(t1_volume_register_dartel_cli.cli)
cli.add_command(t1_volume_tissue_segmentation_cli.cli)


if __name__ == "__main__":
    cli()
