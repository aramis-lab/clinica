from typing import Optional

import click

import clinica.pydra.engine_utils as pydra_utils
import clinica.pydra.t1_linear.t1_linear as pydra_t1_linear
from clinica import option
from clinica.pipelines import cli_param
from clinica.pipelines.engine import clinica_pipeline

pipeline_name = "pydra-t1-linear"


@clinica_pipeline
@click.command(name=pipeline_name, hidden=True)
@cli_param.argument.bids_directory
@cli_param.argument.caps_directory
@option.global_option_group
@option.n_procs
def cli(
    bids_directory: str,
    caps_directory: str,
    n_procs: Optional[int] = None,
) -> None:
    """Affine registration of T1w images to the MNI standard space (Pydra engine).

    Parameters
    ----------
    bids_directory : str
        Path to the input folder containing a BIDS compliant dataset.

    caps_directory : str
        Path to a CAPS folder containing the results of the required pipelines (T1-linear).
        This will also be used for writing the outputs.

    n_procs : int, optional
        Number of processes to use.
    """
    t1_linear_pipeline = pydra_t1_linear.build_core_workflow(
        name="t1-linear-pydra",
        input_dir=bids_directory,
        output_dir=caps_directory,
        parameters={},
    )
    pydra_utils.run(t1_linear_pipeline, n_procs=n_procs)


if __name__ == "__main__":
    cli()
