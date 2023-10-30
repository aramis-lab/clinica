from typing import List, Optional

import click

import clinica.pydra.engine_utils as pydra_utils
import clinica.pydra.t1_volume.tissue_segmentation.pipeline as pydra_t1vol
from clinica import option
from clinica.pipelines import cli_param
from clinica.pipelines.engine import clinica_pipeline

pipeline_name = "pydra-t1-volume-tissue-segmentation"


@clinica_pipeline
@click.command(name=pipeline_name, hidden=True)
@cli_param.argument.bids_directory
@cli_param.argument.caps_directory
@cli_param.option_group.common_pipelines_options
@cli_param.option.subjects_sessions_tsv
@cli_param.option.working_directory
@cli_param.option.yes
@option.global_option_group
@option.n_procs
@cli_param.option_group.advanced_pipeline_options
@cli_param.option.tissue_classes
@cli_param.option.dartel_tissues
@cli_param.option.tissue_probability_maps
@cli_param.option.dont_save_warped_unmodulated
@cli_param.option.save_warped_modulated
def cli(
    bids_directory: str,
    caps_directory: str,
    tissue_classes: List[int] = (1, 2, 3),
    dartel_tissues: List[int] = (1, 2, 3),
    tissue_probability_maps: Optional[str] = None,
    dont_save_warped_unmodulated: bool = False,
    save_warped_modulated: bool = False,
    subjects_sessions_tsv: Optional[str] = None,
    working_directory: Optional[str] = None,
    n_procs: Optional[int] = None,
    yes: bool = False,
) -> None:
    """Affine registration of Flair images to the MNI standard space (Pydra engine)."""

    parameters = {
        "tissue_classes": tissue_classes,
        "dartel_tissues": dartel_tissues,
        "tissue_probability_maps": tissue_probability_maps,
        "save_warped_unmodulated": not dont_save_warped_unmodulated,
        "save_warped_modulated": save_warped_modulated,
        "skip_question": yes,
    }

    t1_volume_tissue_segmentation_pipeline = pydra_t1vol.t1volume_tissue_segmentation(
        name="t1-volume-tissue-segmentation-pydra",
        input_dir=bids_directory,
        output_dir=caps_directory,
        parameters=parameters,
    )
    pydra_utils.run(t1_volume_tissue_segmentation_pipeline, n_procs=n_procs)


if __name__ == "__main__":
    cli()
