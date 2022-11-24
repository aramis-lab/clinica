from typing import List

import click

import clinica.pydra.engine_utils as pydra_utils
from clinica.pipelines import cli_param

pipeline_name = "pydra-t1-volume-register-dartel"


@click.command(name=pipeline_name, hidden=True)
@cli_param.argument.bids_directory
@cli_param.argument.caps_directory
@cli_param.argument.group_label
# @cli_param.option_group.common_pipelines_options
# @cli_param.option.subjects_sessions_tsv
# @cli_param.option.working_directory
# @cli_param.option.n_procs
@cli_param.option_group.advanced_pipeline_options
@cli_param.option.dartel_tissues
def cli(
    bids_directory: str,
    caps_directory: str,
    group_label: str,
    dartel_tissues: List[int] = (1, 2, 3),
) -> None:
    """Inter-subject registration using Dartel (using an existing Dartel template).

    Parameters
    ----------

    group_label : str
        User-defined identifier to target a specific group of subjects. For this pipeline, it is associated to the
        DARTEL template that you had created when running the t1-volume pipeline.

    dartel_tissues : List
        indices for tissues in ["graymatter", 2: "whitematter", 3: "csf", 4: "bone", 5: "softtissue", 6: "background"]
    Returns
    -------
    None

    Notes
    -----
    https://aramislab.paris.inria.fr/clinica/docs/public/latest/Pipelines/T1_Volume/
    """

    import clinica.pydra.t1_volume.register_dartel.pipeline as pydra_t1vol_rd

    parameters = {"group_label": group_label, "dartel_tissues": dartel_tissues}

    t1_volume_tissue_segmentation_pipeline = pydra_t1vol_rd.t1volume_register_dartel(
        name="t1-volume-tissue-segmentation-pydra",
        input_dir=bids_directory,
        output_dir=caps_directory,
        parameters=parameters,
    )
    pydra_utils.run(t1_volume_tissue_segmentation_pipeline)
