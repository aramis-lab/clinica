from typing import List, Optional, Tuple

import click

from clinica import option
from clinica.pipelines import cli_param
from clinica.pipelines.engine import clinica_pipeline

pipeline_name = "pet-volume"


@clinica_pipeline
@click.command(name=pipeline_name)
@cli_param.argument.bids_directory
@cli_param.argument.caps_directory
@cli_param.argument.group_label
@cli_param.argument.acq_label
@cli_param.argument.suvr_reference_region
@cli_param.option_group.pipeline_specific_options
@cli_param.option.pvc_psf_tsv
@cli_param.option_group.common_pipelines_options
@cli_param.option.reconstruction_method
@cli_param.option.subjects_sessions_tsv
@cli_param.option.working_directory
@cli_param.option.yes
@option.global_option_group
@option.n_procs
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
    reconstruction_method: Optional[str] = None,
    pvc_psf_tsv: Optional[str] = None,
    mask_tissues: List[int] = (1, 2, 3),
    mask_threshold: float = 0.3,
    pvc_mask_tissues: List[int] = (1, 2, 3),
    smooth: List[int] = (8, 8),
    subjects_sessions_tsv: Optional[str] = None,
    working_directory: Optional[str] = None,
    n_procs: Optional[int] = None,
    yes: bool = False,
) -> None:
    """SPM-based pre-processing of PET images.

       GROUP_LABEL is an user-defined identifier to target a specific group of subjects.
    For this pipeline, it is associated to the DARTEL template that you had created when running the t1-volume pipeline.

       ACQ_LABEL corresponds to the label given to the PET acquisition, specifying the tracer used.
    Frequently used values are '18FFDG' or '18FAV45'.

       The reference region must be specified to perform intensity normalization.
    Accepted values include: 'pons', 'cerebellumPons', 'pons2', 'cerebellumPons2'.

    Prerequisite: You need to have performed the t1-volume pipeline on your T1-weighted MR images.

    See https://aramislab.paris.inria.fr/clinica/docs/public/latest/Pipelines/PET_Volume/
    """
    from networkx import Graph

    from clinica.pipelines.pet_volume.pet_volume_pipeline import PETVolume
    from clinica.utils.ux import print_end_pipeline

    parameters = {
        "group_label": group_label,
        "acq_label": acq_label,
        "suvr_reference_region": suvr_reference_region,
        "reconstruction_method": reconstruction_method,
        "pvc_psf_tsv": pvc_psf_tsv,
        "mask_tissues": mask_tissues,
        "mask_threshold": mask_threshold,
        "pvc_mask_tissues": pvc_mask_tissues,
        "smooth": smooth,
        "skip_question": yes,
    }

    pipeline = PETVolume(
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


if __name__ == "__main__":
    cli()
