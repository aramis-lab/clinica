from pathlib import Path

import pandas as pd

from clinica.utils.stream import cprint

__all__ = ["find_image_path"]


def find_image_path(
    images: pd.DataFrame,
    source_dir: Path,
    modality: str,
) -> pd.DataFrame:
    """
    For each image, the path to an existing image file or folder is created from image metadata.

    This function adds two columns to the input dataframe: 'Is_Dicom', and 'Path'.

    Args:
        images: List of images metadata
        source_dir: path to the ADNI directory
        modality: Imaging modality

    Returns: Dataframe containing metadata and existing paths
    """
    is_dicom = []
    image_folders = []
    for _, image in images.iterrows():
        image_path, dicom = _find_path_single_image(image, source_dir)
        is_dicom.append(dicom)
        image_folders.append(image_path)
        if image_path == "":
            cprint(
                msg=(
                    f"No {modality} image path found for subject {image.Subject_ID} in visit {image.VISCODE} "
                    f"with image ID {image.Image_ID}"
                ),
                lvl="info",
            )
    images.loc[:, "Is_Dicom"] = pd.Series(is_dicom, index=images.index)
    images.loc[:, "Path"] = pd.Series(image_folders, index=images.index)

    return images


def _find_path_single_image(image: pd.Series, source_dir: Path) -> tuple[str, bool]:
    path_to_sequence = source_dir / str(image["Subject_ID"])
    image_folder_path = ""
    is_dicom = True
    if (
        len(
            (
                found_images := [
                    f
                    for f in path_to_sequence.rglob(f"*_I{image['Image_ID']}.*")
                    if f.is_file()
                ]
            )
        )
        > 0
    ):
        if "dcm" in found_images[0].suffix:
            image_folder_path = str(found_images[0].parent)
        else:
            image_folder_path = str(found_images[0])
            is_dicom = False
    return image_folder_path, is_dicom
