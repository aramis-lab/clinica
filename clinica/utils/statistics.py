"""This module contains utilities for statistics.

Currently, it contains one function to generate TSV file containing mean map based on a parcellation.
"""
from os import PathLike
from pathlib import Path
from typing import Optional, Union

from clinica.utils.atlas import AtlasName, BaseAtlas


def statistics_on_atlas(
    in_normalized_map: Union[str, PathLike],
    atlas: Union[str, AtlasName, BaseAtlas],
    out_file: Optional[Union[str, PathLike]] = None,
) -> str:
    """Compute statistics of a map on an atlas.

    Given an atlas image with a set of ROIs, this function computes the mean of
    a normalized map (e.g. GM segmentation, FA map from DTI, etc.) on each ROI.

    Parameters
    ----------
    in_normalized_map : str
        File containing a scalar image registered on the atlas.

    atlas : BaseAtlas or AtlasName or str
        An atlas with a set of ROI. These ROI are used to compute statistics.
        If a string is given, it is assumed to be the name of the atlas to be used.

    out_file : str, optional
        Name of the output file.

    Returns
    -------
    out_file : str
        TSV file containing the statistics (content of the columns: label,
        mean scalar, std of the scalar', number of voxels).
    """
    import nibabel as nib
    import numpy as np
    import pandas as pd

    from clinica.utils.stream import cprint

    from .atlas import atlas_factory

    atlas = atlas_factory(atlas)
    in_normalized_map = Path(in_normalized_map)
    if not out_file:
        filename, ext = in_normalized_map.stem, in_normalized_map.suffix
        if ext == ".gz":
            filename = Path(filename).stem
        out_file = Path(f"{filename}_statistics_{atlas.name}.tsv").resolve()

    atlas_labels = nib.load(atlas.labels)
    atlas_labels_data = atlas_labels.get_fdata(dtype="float32")

    img = nib.load(in_normalized_map)
    img_data = img.get_fdata(dtype="float32")

    atlas_correspondence = pd.read_csv(atlas.tsv_roi, sep="\t")
    label_name = list(atlas_correspondence.roi_name)
    # TODO create roi_value column in lut_*.txt and remove irrelevant RGB information
    label_value = list(atlas_correspondence.roi_value)

    mean_signal_value = []
    for label in label_value:
        current_mask_label = atlas_labels_data == label
        masked_data = np.array(img_data, copy=True)
        masked_data[np.invert(current_mask_label)] = 0
        mean_signal_value.append(np.sum(masked_data) / np.sum(current_mask_label))

    try:
        data = pd.DataFrame(
            {"label_name": label_name, "mean_scalar": mean_signal_value}
        )
        data.to_csv(out_file, sep="\t", index=True, encoding="utf-8")
    except Exception as e:
        cprint(msg=f"Impossible to save {out_file} with pandas", lvl="error")
        raise e

    return out_file
