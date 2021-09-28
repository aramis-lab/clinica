from typing import Optional

import click

from clinica.pipelines import cli_param

pipeline_name = "compute-atlas"


@click.command(name=pipeline_name)
@cli_param.argument.bids_directory
@cli_param.argument.caps_directory
@cli_param.option.atlas_path
def cli(
    bids_directory: str,
    caps_directory: str,
    atlas_path: Optional[str] = None,
) -> None:
    """Cross-sectional pre-processing of T1w images with FreeSurfer.

    https://aramislab.paris.inria.fr/clinica/docs/public/latest/Pipelines/T1_FreeSurfer/
    """
    from networkx import Graph

    from clinica.utils.ux import print_end_pipeline

    from .compute_atlas_pipeline import ComputeAtlas

    pipeline = ComputeAtlas(
        bids_directory=bids_directory,
        caps_directory=caps_directory,
        atlas_path=atlas_path,
    )

    exec_pipeline = pipeline.run()

    if isinstance(exec_pipeline, Graph):
        print_end_pipeline(
            pipeline_name, pipeline.base_dir, pipeline.base_dir_was_specified
        )


if __name__ == "__main__":
    cli()
