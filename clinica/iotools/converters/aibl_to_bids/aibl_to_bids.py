"""
Convert the AIBL dataset (https://www.aibl.csiro.au/) into BIDS.
"""
from os import PathLike
from typing import List, Optional

from clinica.iotools.abstract_converter import Converter


class AiblToBidsConverter(Converter):
    study_name: str = "AIBL"

    def __init__(
        self,
        source_dataset: PathLike,
        destination_dataset: PathLike,
        clinical_data_directory: PathLike,
        clinical_data_only: bool = False,
        overwrite: bool = False,
    ):
        super().__init__(
            source_dataset,
            destination_dataset,
            clinical_data_directory,
            clinical_data_only,
        )
        self.overwrite = overwrite

    def convert_images(
        self,
        subjects_list_path: Optional[PathLike] = None,
        modalities: Optional[List[str]] = None,
    ) -> None:
        """Conversion of the AIBL imaging data in BIDS.

        overwrite : bool, optional
            Overwrites previously written nifti and json files.
            Default=False.
        """
        from os.path import exists

        from clinica.iotools.converters.aibl_to_bids.utils import (
            Modality,
            paths_to_bids,
        )
        from clinica.utils.stream import cprint

        list_of_created_files = [
            paths_to_bids(
                self.source_dataset,
                self.clinical_data_directory,
                self.destination_dataset,
                modality,
                overwrite=self.overwrite,
            )
            for modality in Modality
        ]
        missing_files = []
        for modality_list in list_of_created_files:
            for file in modality_list:
                if not exists(str(file)):
                    missing_files.append(file)
        if missing_files:
            msg = "The following file were not converted:\n" + "\n".join(missing_files)
            cprint(msg=msg, lvl="warning")

    def convert_clinical_data(
        self, subjects_list_path: Optional[PathLike] = None
    ) -> None:
        """Conversion of the AIBL clinical data in BIDS.

        Parameters
        ----------
        """
        from os.path import join, realpath, split

        from clinica.iotools.converters.aibl_to_bids.utils import (
            create_participants_tsv_file,
            create_scans_tsv_file,
            create_sessions_tsv_file,
        )
        from clinica.utils.stream import cprint

        clinical_spec_path = join(
            split(realpath(__file__))[0], "../../data/clinical_specifications"
        )
        cprint("Creating participants.tsv...")
        create_participants_tsv_file(
            self.destination_dataset,
            clinical_spec_path,
            self.clinical_data_directory,
            delete_non_bids_info=True,
        )
        cprint("Creating sessions files...")
        create_sessions_tsv_file(
            self.destination_dataset, self.clinical_data_directory, clinical_spec_path
        )
        cprint("Creating scans files...")
        create_scans_tsv_file(
            self.destination_dataset, self.clinical_data_directory, clinical_spec_path
        )
        super().convert_clinical_data(subjects_list_path)
