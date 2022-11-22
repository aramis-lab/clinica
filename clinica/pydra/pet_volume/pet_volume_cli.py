from typing import List, Optional, Union

import click

import clinica.pydra.engine_utils as pydra_utils
import clinica.pydra.pet_volume.pipeline as pydra_pet_volume
from clinica.pipelines import cli_param
from clinica.pipelines.cli import cli as run_cli

pipeline_name = "pydra-pet-volume"


@click.command(name=pipeline_name, hidden=True)
@cli_param.argument.bids_directory
@cli_param.argument.caps_directory
@cli_param.argument.group_label
@cli_param.argument.acq_label
@cli_param.argument.suvr_reference_region
@cli_param.option_group.pipeline_specific_options
@cli_param.option.pvc_psf_tsv
@cli_param.option_group.common_pipelines_options
@cli_param.option_group.advanced_pipeline_options
@cli_param.option_group.option(
    "-mask",
    "--mask_tissues",
    multiple=True,
    type=click.IntRange(1, 7),
    default=[1, 2, 3],
    help=(
        "Tissue classes (1: gray matter (GM), 2: white matter (WM), 3: cerebrospinal "
        "fluid (CSF), 4: bone, 5: soft-tissue, 6: background) to use "
        "for masking the PET image (default: GM, WM and CSF i.e. --mask_tissues 1 2 3)."
    ),
)
@cli_param.option_group.option(
    "-threshold",
    "--mask_threshold",
    default=0.3,
    show_default=True,
    help="Value used as threshold to binarize the tissue maps.",
)
@cli_param.option_group.option(
    "-pvc_mask",
    "--pvc_mask_tissues",
    multiple=True,
    type=click.IntRange(1, 7),
    default=[1, 2, 3],
    help=(
        "Tissue classes (1: gray matter (GM), 2: white matter (WM), 3: cerebrospinal "
        "fluid (CSF), 4: bone, 5: soft-tissue, 6: background) to use "
        "as mask for PVC (default: GM, WM and CSF i.e. --pvc_mask_tissues 1 2 3)."
    ),
)
@cli_param.option.smooth
def cli(
    bids_directory: str,
    caps_directory: str,
    group_label: str,
    acq_label: str,
    suvr_reference_region: Optional[str] = None,
    pvc_psf_tsv: Optional[str] = None,
    mask_tissues: List[int] = (1, 2, 3),
    mask_threshold: float = 0.3,
    pvc_mask_tissues: List[int] = (1, 2, 3),
    smooth: Union[float, List[float], List[List[float]]] = 8.0,
) -> None:
    """SPM-based pre-processing of PET images (Pydra engine).

    GROUP_LABEL is an user-defined identifier to target a specific group of subjects.

    This pipeline is associated to the DARTEL template that you had created when
    running the t1-volume pipeline.

    ACQ_LABEL corresponds to the label given to the PET acquisition, specifying the
    tracer used.
    Frequently used values are '18FFDG' or '18FAV45'.

    The reference region must be specified to perform intensity normalization.
    Accepted values include: 'pons', 'cerebellumPons', 'pons2', 'cerebellumPons2'.

    Prerequisite: You need to have performed the t1-volume pipeline on your T1-weighted MR images.

    See https://aramislab.paris.inria.fr/clinica/docs/public/latest/Pipelines/PET_Volume/
    """
    parameters = {
        "group_label": group_label,
        "acq_label": acq_label,
        "suvr_reference_region": suvr_reference_region,
        "pvc_psf_tsv": pvc_psf_tsv,
        "mask_tissues": mask_tissues,
        "mask_threshold": mask_threshold,
        "pvc_mask_tissues": pvc_mask_tissues,
        "smooth": smooth,
    }
    pipeline = pydra_pet_volume.build_core_workflow(
        name="pet-volume-pydra",
        input_dir=bids_directory,
        output_dir=caps_directory,
        parameters=parameters,
    )
    pydra_utils.run(pipeline)


run_cli.add_command(cli)

if __name__ == "__main__":
    cli()
