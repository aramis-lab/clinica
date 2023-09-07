"""Convert the GENFI dataset into BIDS."""

from os import PathLike
from pathlib import Path
from typing import List, Optional

from clinica.iotools.abstract_converter import Converter


class GenfiToBidsConverter(Converter):
    """BIDS converter for Genfi dataset.

    Attributes
    ----------
    gif: bool
        If True, indicates the user wants to have the values of the gif parcellation
    """

    study_name: str = "GENFI"

    def __init__(
        self,
        source_dataset: PathLike,
        destination_dataset: PathLike,
        clinical_data_directory: PathLike,
        clinical_data_only: bool = False,
        gif: bool = False,
        clinical_data: bool = True,
    ):
        super().__init__(
            source_dataset,
            destination_dataset,
            clinical_data_directory,
            clinical_data_only,
        )
        self.gif = gif
        self.clinical_data = clinical_data

    @property
    def readme_path(self) -> Path:
        import os

        return Path(
            os.path.join(
                os.path.dirname(
                    os.path.dirname(
                        os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
                    )
                ),
                "docs",
                "Converters",
                "GENFItoBIDS.md",
            )
        )

    def convert_images(
        self,
        subjects_list_path: Optional[PathLike] = None,
        modalities: Optional[List[str]] = None,
    ) -> None:
        """Convert the entire dataset to BIDS.

        Scans available files in the path_to_dataset,
        identifies the patients that have images described by the JSON file,
        converts the image with the highest quality for each category.
        """
        from .genfi_to_bids_utils import (
            complete_clinical_data,
            dataset_to_bids,
            find_clinical_data,
            intersect_data,
            merge_imaging_data,
            read_imaging_data,
            write_bids,
        )

        # read the clinical data files
        if self.clinical_data:
            df_demographics, df_imaging, df_clinical = find_clinical_data(
                self.clinical_data_directory
            )
        # makes a df of the imaging data
        imaging_data = read_imaging_data(self.source_dataset)

        # complete the data extracted
        imaging_data = merge_imaging_data(imaging_data)
        # complete clinical data
        if self.clinical_data:
            df_clinical_complete = complete_clinical_data(
                df_demographics, df_imaging, df_clinical
            )

        # intersect the data
        if self.clinical_data:
            df_complete = intersect_data(imaging_data, df_clinical_complete)
        else:
            df_complete = imaging_data
        # build the tsv
        results = dataset_to_bids(df_complete, self.gif)
        write_bids(
            to=self.destination_dataset,
            participants=results["participants"],
            sessions=results["sessions"],
            scans=results["scans"],
        )
