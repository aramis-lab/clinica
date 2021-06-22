from typing import Optional

import click

from clinica.pipelines import cli_param

pipeline_name = "statistics-volume-correction"


@click.command(name=pipeline_name)
@cli_param.argument.caps_directory
@click.argument("t_map", type=click.Path(exists=True, resolve_path=True))
@click.argument("height_threshold", type=float)
@click.argument("FWEp", type=float)
@click.argument("FDRp", type=float)
@click.argument("FWEc", type=float)
@click.argument("FDRc", type=float)
@cli_param.option_group.standard_options
@cli_param.option.working_directory
@cli_param.option.n_procs
@click.option(
    "-nc",
    "--n_cuts",
    default=8,
    show_default=True,
    help="Number of cuts along each direction",
)
def cli(
    caps_directory: str,
    t_map: str,
    height_threshold: str,
    fwep: float,
    fdrp: float,
    fwec: float,
    fdrc: float,
    n_cuts: int = 8,
    working_directory: Optional[str] = None,
    n_procs: Optional[int] = None,
) -> None:
    """Statistical correction of statistics-volume pipeline:

    https://aramislab.paris.inria.fr/clinica/docs/public/latest/Pipelines/Stats_Volume/
    """
    from networkx import Graph

    from clinica.utils.ux import print_end_pipeline

    from .statistics_volume_correction_pipeline import StatisticsVolumeCorrection

    parameters = {
        "t_map": t_map,
        "height_threshold": height_threshold,
        "FWEp": fwep,
        "FDRp": fdrp,
        "FWEc": fwec,
        "FDRc": fdrc,
        "n_cuts": n_cuts,
    }

    pipeline = StatisticsVolumeCorrection(
        caps_directory=caps_directory,
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


if __name__ == "__main__":
    cli()
