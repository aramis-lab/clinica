import click

import clinica.pydra.dwi_preprocessing_using_t1 as pydra_dwi_preprocessing_using_t1_pipeline
import clinica.pydra.engine_utils as pydra_utils
from clinica.pipelines import cli_param

pipeline_name = "pydra-dwi-preprocessing-using-t1"


@click.command(name=pipeline_name, hidden=True)
@cli_param.argument.bids_directory
@cli_param.argument.caps_directory
def cli(
    bids_directory: str,
    caps_directory: str,
) -> None:
    """Affine registration of Flair images to the MNI standard space (Pydra engine)."""

    dwi_preprocessing_using_t1_pipeline = (
        pydra_dwi_preprocessing_using_t1_pipeline.build_core_workflow(
            name="pydra-dwi-preprocessing-using-t1",
            input_dir=bids_directory,
            output_dir=caps_directory,
        )
    )
    pydra_utils.run(dwi_preprocessing_using_t1_pipeline)


if __name__ == "__main__":
    cli()
