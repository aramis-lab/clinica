import abc
from os import PathLike
from pathlib import Path
from typing import Dict, List, Optional


class Converter:
    """Converter base class for conversion of datasets to BIDS.

    Attributes
    ----------
    study_name : str
        The name of the study (ex: UKB, OASIS, ADNI...).

    source_dataset : Path
        The path to the input raw dataset to be converted.

    destination_dataset : Path
        The path to the destination folder in which the BIDS
        dataset converted from the source dataset should be
        written. If the folders do not exist, they will be
        created.

    clinical_data_directory : Path
        The path to the input folder with the clinical data.

    clinical_data_only : bool
        If set to True, the converter will not try to convert images
        which requires resources to be performed, only the clinical
        data will be converted then.
        If set to False, the converter will try to convert everything.
    """

    study_name: str

    def __init__(
        self,
        source_dataset: PathLike,
        destination_dataset: PathLike,
        clinical_data_directory: PathLike,
        clinical_data_only: bool = False,
    ):
        self.source_dataset = Path(source_dataset)
        self.destination_dataset = Path(destination_dataset)
        if not self.destination_dataset.exists():
            self.destination_dataset.mkdir(parents=True, exist_ok=True)
        self.clinical_data_directory = Path(clinical_data_directory)
        self.clinical_data_only = clinical_data_only

    @property
    def link(self) -> str:
        """The URL to download the raw source dataset which is the input of the converter."""
        from clinica.iotools.bids_utils import parse_url

        try:
            return parse_url(self.readme_path)[0]
        except IndexError:
            raise ValueError("Could not parse URL of dataset.")

    @property
    def description(self) -> str:
        """Short description of the dataset.

        For more information please refer to the online documentation.
        """
        from clinica.iotools.bids_utils import parse_description

        start_line, end_line = 4, 5
        try:
            return parse_description(self.readme_path, start_line, end_line)
        except IndexError:
            raise ValueError("Could not parse description for dataset.")

    @property
    def readme_path(self) -> Path:
        """The path to the corresponding README file in the docs.

        This is used internally by the converter to parse the link and description.
        """
        from clinica.utils.filemanip import get_parent

        return (
            get_parent(__file__, 3)
            / "docs"
            / "Converters"
            / f"{self.study_name}toBIDS.md"
        )

    @property
    def readme_data(self) -> Dict[str, str]:
        return {"link": self.link, "desc": self.description}

    def check_dependencies(self) -> None:
        """Check the dependencies required to run the converter."""
        from clinica.utils.check_dependency import check_dcm2niix

        check_dcm2niix()

    def convert(
        self,
        subjects_list_path: Optional[PathLike] = None,
        modalities: Optional[List[str]] = None,
    ) -> None:
        self.check_dependencies()
        if not self.clinical_data_only:
            self.convert_images(subjects_list_path, modalities)
        self.convert_clinical_data(subjects_list_path)

    @abc.abstractmethod
    def convert_images(
        self,
        subjects_list_path: Optional[PathLike] = None,
        modalities: Optional[List[str]] = None,
    ) -> None:
        """Convert the image data.

        Parameters
        ----------
        subjects_list_path : str, optional
            If specified, restrict the processing to the subjects specified
            in the corresponding file.
            If not specified, all subjects will be handled.
        """
        pass

    def convert_clinical_data(
        self, subjects_list_path: Optional[PathLike] = None
    ) -> None:
        """Convert the clinical data.

        Parameters
        ----------
        subjects_list_path : str, optional
            If specified, restrict the processing to the subjects specified
            in the corresponding file.
            If not specified, all subjects will be handled.
        """
        self.write_modality_agnostic_files()

    def write_modality_agnostic_files(self) -> None:
        from clinica.iotools.bids_utils import write_modality_agnostic_files
        from clinica.utils.stream import cprint

        cprint("Creating modality agnostic files...")
        write_modality_agnostic_files(
            study_name=self.study_name,
            readme_data=self.readme_data,
            bids_dir=self.destination_dataset,
        )
