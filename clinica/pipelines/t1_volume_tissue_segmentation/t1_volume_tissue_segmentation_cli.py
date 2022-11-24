from typing import List, Optional

import click

from clinica.pipelines import cli_param
from clinica.pipelines.cli import cli as run_cli

pipeline_name = "t1-volume-tissue-segmentation"


@click.command(name=pipeline_name)
@cli_param.argument.bids_directory
@cli_param.argument.caps_directory
@cli_param.option_group.common_pipelines_options
@cli_param.option.subjects_sessions_tsv
@cli_param.option.working_directory
@cli_param.option.n_procs
@cli_param.option.yes
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
    """Tissue segmentation, bias correction and spatial normalization to MNI space of T1w images with SPM.

    https://aramislab.paris.inria.fr/clinica/docs/public/latest/Pipelines/T1_Volume/
    """
    from networkx import Graph

    from clinica.utils.ux import print_end_pipeline

    from .t1_volume_tissue_segmentation_pipeline import T1VolumeTissueSegmentation

    parameters = {
        "tissue_classes": tissue_classes,
        "dartel_tissues": dartel_tissues,
        "tissue_probability_maps": tissue_probability_maps,
        "save_warped_unmodulated": not dont_save_warped_unmodulated,
        "save_warped_modulated": save_warped_modulated,
        "skip_question": yes,
    }

    pipeline = T1VolumeTissueSegmentation(
        bids_directory=bids_directory,
        caps_directory=caps_directory,
        tsv_file=subjects_sessions_tsv,
        base_dir=working_directory,
        parameters=parameters,
        name=pipeline_name,
    )

    exec_pipeline = (
        pipeline.run(plugin="MultiProc", plugin_args={"n_procs": n_procs})
        if n_procs
        else pipeline.run()
    )

    if isinstance(exec_pipeline, Graph):
        print_end_pipeline(
            pipeline_name, pipeline.base_dir, pipeline.base_dir_was_specified
        )


run_cli.add_command(cli)

if __name__ == "__main__":
    cli()
