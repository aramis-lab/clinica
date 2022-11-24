from typing import List, Optional, Tuple

import click

import clinica.pydra.engine_utils as pydra_utils
import clinica.pydra.t1_volume.dartel2mni.pipeline as pydra_t1_vol_dartel2mni
from clinica.pipelines import cli_param
from clinica.pipelines.engine import clinica_pipeline

pipeline_name = "pydra-t1-volume-dartel2mni"


@clinica_pipeline
@click.command(name=pipeline_name, hidden=True)
@cli_param.argument.bids_directory
@cli_param.argument.caps_directory
def cli(
    bids_directory: str,
    caps_directory: str,
    group_label: str,
    smooth: List[int] = (8,),
    tissues: List[int] = (1, 2, 3),
    modulate: bool = True,
    voxel_size: Tuple[float, float, float] = (1.5, 1.5, 1.5),
    subjects_sessions_tsv: Optional[str] = None,
    working_directory: Optional[str] = None,
    n_procs: Optional[int] = None,
) -> None:
    """
    Register DARTEL template to MNI space (Pydra engine).

    GROUP_LABEL is an user-defined identifier to target a specific group of subjects.
    For this pipeline, it is associated to the DARTEL template that you had created
    when running the t1-volume pipeline.
    https://aramislab.paris.inria.fr/clinica/docs/public/latest/Pipelines/T1_Volume/
    """
    parameters = {
        "group_label": group_label,
        "tissues": tissues,
        "voxel_size": voxel_size,
        "modulate": modulate,
        "smooth": smooth,
    }
    pipeline = pydra_t1_vol_dartel2mni.build_core_workflow(
        name=pipeline_name,
        input_dir=bids_directory,
        output_dir=caps_directory,
        parameters=parameters,
    )
    pydra_utils.run(pipeline)


if __name__ == "__main__":
    cli()
