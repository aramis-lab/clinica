"""Convert the UKB dataset into BIDS."""

from os import PathLike
from typing import List, Optional

from clinica.iotools.abstract_converter import Converter


class UkbToBidsConverter(Converter):
    study_name: str = "UKB"
    # link: str = "https://www.ukbiobank.ac.uk/",
    # description: str = (
    #    "UK Biobank is a large-scale biomedical database and research resource, "
    #    "containing in-depth genetic and health information from half a million UK "
    #    "participants. The database is regularly augmented with additional data and "
    #    "is globally accessible to approved researchers undertaking vital research "
    #    "into the most common and life-threatening diseases. It is a major contributor "
    #    "to the advancement of modern medicine it and has led to the discovery of "
    #    "several scientific advances and numerous treatments to improve human health."
    # )

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
        from .ukb_utils import (
            complete_clinical,
            dataset_to_bids,
            find_clinical_data,
            intersect_data,
            read_imaging_data,
            write_bids,
        )

        # read the clinical data files
        df_clinical = find_clinical_data(self.clinical_data_directory)

        # makes a df of the imaging data
        imaging_data = read_imaging_data(self.source_dataset)

        # intersect the data
        df_clinical = intersect_data(imaging_data, df_clinical)

        # complete clinical data
        df_clinical = complete_clinical(df_clinical)

        # build the tsv
        result = dataset_to_bids(df_clinical)

        write_bids(
            to=self.destination_dataset,
            participants=result["participants"],
            sessions=result["sessions"],
            scans=result["scans"],
            dataset_directory=self.source_dataset,
        )
