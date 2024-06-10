from typing import Optional, Sequence

import numpy as np
import pandas as pd

__all__ = ["select_image_qc"]


def select_image_qc(id_list: Sequence[str], mri_qc_subj: pd.DataFrame) -> Optional[int]:
    """Select image from several scans according to QC.

    Args:
        id_list: List of images identifiers to choose among
        mri_qc_subj: Dataframe containing list of QC of scans for the subject

    Returns: int - Chosen image identifier
    """
    if len(id_list) == 0:
        return None

    image_ids = ["I" + str(imageuid) for imageuid in id_list]
    int_ids = [int(imageuid) for imageuid in id_list]
    images_qc = mri_qc_subj[mri_qc_subj.loni_image.isin(image_ids)]

    if images_qc.empty:
        return max(int_ids)

    if np.sum(images_qc.series_selected) == 1:
        selected_image = (
            images_qc[images_qc.series_selected == 1].iloc[0].loni_image[1:]
        )
    else:
        images_not_rejected = images_qc[images_qc.series_quality < 4]

        if images_not_rejected.empty:
            # There are no images that passed the qc,
            # so we'll try to see if there are other images without qc.
            # Otherwise, return None.
            qc_ids = set([int(qc_id[1:]) for qc_id in images_qc.loni_image.unique()])
            no_qc_ids = list(set(int_ids) - qc_ids)

            if len(no_qc_ids) == 0:
                return None
            else:
                return max(no_qc_ids)

        # We select the image with the best (lower) QC.
        # If no positive QC available we choose image with -1 (not performed)
        series_quality = [
            q if q > 0 else 4 for q in list(images_not_rejected.series_quality)
        ]
        best_q = np.amin(series_quality)
        if best_q == 4:
            best_q = -1
        images_best_qc = images_not_rejected[
            images_not_rejected.series_quality == best_q
        ]
        if len(images_best_qc) == 1:
            selected_image = images_best_qc.iloc[0].loni_image[1:]
        else:
            best_ids = [int(x[1:]) for x in images_best_qc.loni_image.unique()]
            selected_image = max(best_ids)

    return int(selected_image)
