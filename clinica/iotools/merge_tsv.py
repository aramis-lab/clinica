from os import PathLike
from typing import Iterable, Optional, Union

__all__ = ["merge_tsv"]


def merge_tsv(
    bids_directory: Union[str, PathLike],
    output_tsv: Union[str, PathLike],
    caps_directory: Optional[Union[str, PathLike]] = None,
    pipelines: Optional[Iterable[str]] = None,
    volume_atlas_selection: Optional[Iterable[str]] = None,
    freesurfer_atlas_selection: Optional[Iterable[str]] = None,
    pvc_restriction: Optional[int] = None,
    pet_tracers_selection: Optional[Iterable[str]] = None,
    group_selection: Optional[Iterable[str]] = None,
    subjects_sessions_tsv: Optional[Union[str, PathLike]] = None,
    ignore_scan_files: bool = False,
    ignore_session_scan_files: bool = False,
):
    """Merge clinical data into a single TSV file.

    Parameters
    ----------
    bids_directory : PathLike
        The path to the BIDS directory.

    output_tsv : PathLike
        The path to the output TSV file.

    caps_directory : PathLike, optional

    pipelines : Iterable of str, optional
        The names of the pipelines to consider.
        By default, all pipelines will be considered.

    volume_atlas_selection : Iterable of str, optional
        The volume atlases to consider.
        By default, all available atlases will be considered.

    freesurfer_atlas_selection : Iterable of str, optional
        The Freesurfer atlases to consider.
        By default, all available atlases will be considered.

    pvc_restriction : int, optional

    pet_tracers_selection : Iterable of str, optional
        The PET tracers to consider.

    group_selection : Iterable of str, optional
        The group labels to consider.
        By default, all available groups will be considered.

    subjects_sessions_tsv : PathLike, optional

    ignore_scan_files : bool, optional

    ignore_session_scan_files : bool, optional
    """
    from clinica.utils.inputs import check_bids_folder

    from .data_handling import create_merge_file

    check_bids_folder(bids_directory)

    atlas_selection = []
    if volume_atlas_selection is not None:
        atlas_selection += volume_atlas_selection
    if freesurfer_atlas_selection is not None:
        atlas_selection += freesurfer_atlas_selection
    if group_selection == ():
        group_selection = None
    if pet_tracers_selection == ():
        pet_tracers_selection = None

    create_merge_file(
        bids_directory,
        output_tsv,
        caps_dir=caps_directory,
        pipelines=pipelines,
        ignore_scan_files=ignore_scan_files,
        ignore_sessions_files=ignore_session_scan_files,
        atlas_selection=atlas_selection,
        pvc_restriction=pvc_restriction,
        tsv_file=subjects_sessions_tsv,
        group_selection=group_selection,
        tracers_selection=pet_tracers_selection,
    )
