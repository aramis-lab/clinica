import click

import clinica.pydra.engine_utils as pydra_utils
import clinica.pydra.t1_volume.tissue_segmentation.pipeline as pydra_t1vol
from clinica.pipelines import cli_param

pipeline_name = "pydra-t1vol-ts"


@click.command(name=pipeline_name, hidden=True)
@cli_param.argument.bids_directory
@cli_param.argument.caps_directory
def cli(
    bids_directory: str,
    caps_directory: str,
) -> None:
    """Affine registration of Flair images to the MNI standard space (Pydra engine)."""

    t1_volume_tissue_segmentation_pipeline = pydra_t1vol.t1volume_tissue_segmentation(
        name="t1-volume-tissue-segmentation-pydra",
        input_dir=bids_directory,
        output_dir=caps_directory,
    )
    pydra_utils.run(t1_volume_tissue_segmentation_pipeline)


if __name__ == "__main__":
    cli()
