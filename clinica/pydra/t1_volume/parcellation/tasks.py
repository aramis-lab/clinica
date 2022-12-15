from pydra.mark import annotate, task


@task
@annotate({"return": {"atlas_statistics_path": list}})
def atlas_statistics_task(in_image: str, atlas_list: list) -> list:
    """Pydra task for Generate regional measure from atlas_list in TSV files.

    .. note::
        Please refer to the documentation of function
        `clinica.pipelines.t1_volume_parcellation.t1_volume_parcellation_utils.atlas_statistics`.
    """
    from clinica.pipelines.t1_volume_parcellation.t1_volume_parcellation_utils import (
        atlas_statistics,
    )

    return atlas_statistics(in_image, atlas_list)
