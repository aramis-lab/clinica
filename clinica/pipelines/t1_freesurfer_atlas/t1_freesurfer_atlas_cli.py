from typing import Optional

import click

from clinica import option
from clinica.pipelines import cli_param

pipeline_name = "compute-atlas"


@click.command(name=pipeline_name)
@cli_param.argument.caps_directory
@cli_param.option.atlas_path
@option.global_option_group
@option.n_procs
def cli(
    caps_directory: str,
    atlas_path: Optional[str] = None,
    n_procs: Optional[int] = None,
) -> None:
    """Projection of the results of t1-freesurfer on another atlas.

    https://aramislab.paris.inria.fr/clinica/docs/public/latest/Pipelines/T1_FreeSurfer/
    """
    from networkx import Graph

    from clinica.utils.ux import print_end_pipeline

    from .t1_freeesurfer_atlas_pipeline import T1FreeSurferAtlas

    pipeline = T1FreeSurferAtlas(
        caps_directory=caps_directory,
        atlas_path=atlas_path,
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


if __name__ == "__main__":
    cli()
