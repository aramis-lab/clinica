from dataclasses import dataclass
from os import PathLike
from pathlib import Path
from typing import List, Union

import nibabel as nib
import numpy as np

from clinica.utils.image import HemiSphere

__all__ = [
    "FreeSurferAnnotationImageSingleHemisphere",
    "FreeSurferAnnotationImage",
    "FreeSurferAnnotation",
]


@dataclass
class FreeSurferAnnotationImageSingleHemisphere:
    path: Path

    @classmethod
    def from_raw(cls, path: Union[str, PathLike]):
        path = Path(path)
        if path.suffix != ".annot":
            raise ValueError(
                f"File {path} does not define a proper FreeSurfer annotation. "
                "It is expected to be in '.annot' file format."
            )
        return cls(path)


@dataclass
class FreeSurferAnnotationImage:
    left: FreeSurferAnnotationImageSingleHemisphere
    right: FreeSurferAnnotationImageSingleHemisphere

    @classmethod
    def from_raw(cls, left: Union[str, PathLike], right: Union[str, PathLike]):
        return cls(
            FreeSurferAnnotationImageSingleHemisphere.from_raw(left),
            FreeSurferAnnotationImageSingleHemisphere.from_raw(right),
        )


@dataclass
class FreeSurferAnnotation:
    region_names: List[str]
    left_annotations: np.ndarray
    right_annotations: np.ndarray
    left_ctab: np.ndarray
    right_ctab: np.ndarray

    @classmethod
    def from_annotation_image(cls, image: FreeSurferAnnotationImage):
        left_annotations, left_ctab, region_names = nib.freesurfer.io.read_annot(
            image.left.path, orig_ids=False
        )
        right_annotations, right_ctab, _ = nib.freesurfer.io.read_annot(
            image.right.path, orig_ids=False
        )
        return cls(
            [r.astype(str) for r in region_names],
            left_annotations,
            right_annotations,
            left_ctab,
            right_ctab,
        )

    @property
    def n_regions(self) -> int:
        return len(self.region_names)

    def get_annotation(self, hemisphere: Union[str, HemiSphere]) -> np.ndarray:
        hemisphere = HemiSphere(hemisphere)
        if hemisphere == HemiSphere.LEFT:
            return self.left_annotations
        return self.right_annotations

    def get_number_of_vertices(self, hemisphere: Union[str, HemiSphere]) -> int:
        hemisphere = HemiSphere(hemisphere)
        if hemisphere == HemiSphere.LEFT:
            return self.left_annotations.shape[0]
        return self.right_annotations.shape[0]

    def replace_minus_one_annotation_with_zero(self):
        self.left_annotations[self.left_annotations == -1] = 0
        self.right_annotations[self.right_annotations == -1] = 0

    def get_lateralized_region_names(self, left_first: bool = True) -> List[str]:
        hemispheres = (
            (HemiSphere.LEFT, HemiSphere.RIGHT)
            if left_first
            else (HemiSphere.RIGHT, HemiSphere.LEFT)
        )
        return [
            f"{region_name}_{hemisphere.value}"
            for region_name in self.region_names
            for hemisphere in hemispheres
        ]
