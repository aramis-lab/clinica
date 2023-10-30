from typing import List, Optional, Tuple

import click

from clinica import option
from clinica.pipelines import cli_param
from clinica.pipelines.engine import clinica_pipeline

pipeline_name = "t1-volume-existing-template"


@clinica_pipeline
@click.command(name=pipeline_name)
@cli_param.argument.bids_directory
@cli_param.argument.caps_directory
@cli_param.argument.group_label
@cli_param.option_group.pipeline_specific_options
@cli_param.option.smooth
@cli_param.option_group.common_pipelines_options
@cli_param.option.subjects_sessions_tsv
@cli_param.option.working_directory
@option.global_option_group
@option.n_procs
@cli_param.option_group.advanced_pipeline_options
@cli_param.option.tissue_classes
@cli_param.option.tissue_probability_maps
@cli_param.option.dont_save_warped_unmodulated
@cli_param.option.save_warped_modulated
@cli_param.option.dartel_tissues
@cli_param.option.tissues
@cli_param.option.modulate
@cli_param.option.voxel_size
@click.pass_context
def cli(
    ctx: click.Context,
    bids_directory: str,
    caps_directory: str,
    group_label: str,
    smooth: List[int] = (8,),
    tissue_classes: List[int] = (1, 2, 3),
    tissue_probability_maps: Optional[str] = None,
    dont_save_warped_unmodulated: bool = False,
    save_warped_modulated: bool = False,
    dartel_tissues: List[int] = (1, 2, 3),
    tissues: List[int] = (1, 2, 3),
    modulate: bool = True,
    voxel_size: Tuple[float, float, float] = (1.5, 1.5, 1.5),
    subjects_sessions_tsv: Optional[str] = None,
    working_directory: Optional[str] = None,
    n_procs: Optional[int] = None,
) -> None:
    """Volume-based processing of T1-weighted MR images using an existing DARTEL template.

       GROUP_LABEL is an user-defined identifier to target a specific group of subjects.

    https://aramislab.paris.inria.fr/clinica/docs/public/latest/Pipelines/T1_Volume/
    """
    from clinica.utils.stream import cprint

    from ..t1_volume_dartel2mni import t1_volume_dartel2mni_cli
    from ..t1_volume_parcellation import t1_volume_parcellation_cli
    from ..t1_volume_register_dartel import t1_volume_register_dartel_cli
    from ..t1_volume_tissue_segmentation import t1_volume_tissue_segmentation_cli

    cprint(
        "The t1-volume-existing-template pipeline is divided into 4 parts:\n"
        "\tt1-volume-tissue-segmentation pipeline: "
        "Tissue segmentation, bias correction and spatial normalization to MNI space\n"
        "\tt1-volume-register-dartel pipeline: "
        "Inter-subject registration using an existing DARTEL template\n"
        "\tt1-volume-dartel2mni pipeline: "
        "DARTEL template to MNI\n"
        "\tt1-volume-parcellation pipeline: "
        "Atlas statistics"
    )

    cprint("Part 1/4: Running t1-volume-segmentation pipeline")

    ctx.invoke(
        t1_volume_tissue_segmentation_cli.cli,
        bids_directory=bids_directory,
        caps_directory=caps_directory,
        tissue_classes=tissue_classes,
        dartel_tissues=dartel_tissues,
        tissue_probability_maps=tissue_probability_maps,
        dont_save_warped_unmodulated=dont_save_warped_unmodulated,
        save_warped_modulated=save_warped_modulated,
        subjects_sessions_tsv=subjects_sessions_tsv,
        working_directory=working_directory,
        n_procs=n_procs,
    )

    cprint("Part 2/4: Running t1-volume-register-dartel pipeline")
    ctx.invoke(
        t1_volume_register_dartel_cli.cli,
        bids_directory=bids_directory,
        caps_directory=caps_directory,
        group_label=group_label,
        tissues=tissues,
        subjects_sessions_tsv=subjects_sessions_tsv,
        working_directory=working_directory,
        n_procs=n_procs,
    )

    cprint("Part 3/4: Running t1-volume-dartel2mni pipeline")
    ctx.invoke(
        t1_volume_dartel2mni_cli.cli,
        bids_directory=bids_directory,
        caps_directory=caps_directory,
        group_label=group_label,
        smooth=smooth,
        tissues=tissues,
        modulate=modulate,
        voxel_size=voxel_size,
        subjects_sessions_tsv=subjects_sessions_tsv,
        working_directory=working_directory,
        n_procs=n_procs,
    )

    cprint("Part 4/4: Running t1-volume-parcellation pipeline")
    ctx.invoke(
        t1_volume_parcellation_cli.cli,
        caps_directory=caps_directory,
        group_label=group_label,
        subjects_sessions_tsv=subjects_sessions_tsv,
        working_directory=working_directory,
        n_procs=n_procs,
    )


if __name__ == "__main__":
    pass
