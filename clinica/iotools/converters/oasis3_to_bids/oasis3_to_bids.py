"""Convert the NIFD dataset into BIDS."""

from os import PathLike
from typing import List, Optional

from clinica.iotools.abstract_converter import Converter


class Oasis3ToBidsConverter(Converter):
    study_name: str = "OASIS-3"
    link: str = "https://www.oasis-brains.org/#access"
    description: str = (
        (
            "OASIS-3 is a retrospective compilation of data for 1378 participants that were "
            "collected across several ongoing projects through the WUSTL Knight ADRC over the "
            "course of 30years. Participants include 755 cognitively normal adults and 622 "
            "individuals at various stages of cognitive decline ranging in age from 42-95yrs. "
            "All participants were assigned a new random identifier and all dates were removed and "
            "normalized to reflect days from entry into study. The dataset contains 2842 MR "
            "sessions which include T1w, T2w, FLAIR, ASL, SWI, time of flight, resting-state BOLD, "
            "and DTI sequences. Many of the MR sessions are accompanied by volumetric segmentation "
            "files produced through FreeSurfer processing. PET imaging from different tracers, PIB, "
            "AV45, and FDG, totaling over 2157 raw imaging scans and the accompanying post-processed "
            "files from the Pet Unified Pipeline (PUP) are also available in OASIS-3."
        ),
    )

    def convert_images(
        self,
        subjects_list_path: Optional[PathLike] = None,
        modalities: Optional[List[str]] = None,
    ) -> None:
        """Convert the entire dataset in BIDS.

        Scans available files in the path_to_dataset,
        identifies the patients that have images described by the JSON file,
        converts the image with the highest quality for each category.
        """
        from .oasis3_utils import (
            dataset_to_bids,
            intersect_data,
            read_clinical_data,
            read_imaging_data,
            write_bids,
        )

        # read the clinical data files
        dict_df = read_clinical_data(self.clinical_data_directory)

        # makes a df of the imaging data
        imaging_data = read_imaging_data(self.source_dataset)

        # intersect the data
        imaging_data, df_small = intersect_data(imaging_data, dict_df)

        # build the tsv
        participants, sessions, scans = dataset_to_bids(imaging_data, df_small)

        write_bids(
            to=self.destination_dataset,
            participants=participants,
            sessions=sessions,
            scans=scans,
            dataset_directory=self.source_dataset,
        )
