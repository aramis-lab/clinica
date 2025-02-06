"""This module handles Clinica exceptions."""


class ClinicaException(Exception):
    """Base class for Clinica exceptions."""


class ClinicaMissingDependencyError(ClinicaException):
    """Base class for Clinica dependencies errors."""


class ClinicaEnvironmentVariableError(ClinicaException):
    """Something is wrong with an environment variable managed by Clinica."""


class ClinicaDatasetError(ClinicaException):
    """Base class for errors related to Datasets."""


class ClinicaBIDSError(ClinicaDatasetError):
    """Base class for BIDS dataset errors."""


class ClinicaCAPSError(ClinicaDatasetError):
    """Base class for CAPS dataset errors."""


class ClinicaExistingDatasetError(ClinicaDatasetError):
    """Base class for existing dataset errors."""

    def __init__(self, dataset_folder_path):
        super().__init__(
            f"Dataset located at {dataset_folder_path} already contain some files."
        )


class ClinicaImageDimensionError(ClinicaException):
    """Base class for errors linked to image dimensions."""


class ClinicaParserError(ClinicaException):
    """Base class for parser errors."""


class ClinicaXMLParserError(ClinicaParserError):
    """Base class for XML parser errors."""


class ClinicaPipelineConfigurationError(ClinicaException):
    """Base class for configuration errors of clinica pipelines."""


class ClinicaInconsistentDatasetError(ClinicaException):
    """Base class for inconsistent datasets errors."""

    def __init__(self, cross_subj: list):
        max_subjects_displayed = 50
        bound = min(len(cross_subj), max_subjects_displayed)
        self.msg = (
            f"{len(cross_subj)} subjects in your BIDS folder did not respect the longitudinal organisation "
            f"from BIDS specification.\nThe subjects concerned are (showing only the first {bound}):\n\t- "
        )
        self.msg += "\n\t- ".join(cross_subj[:bound])
        self.msg += (
            f"\nClinica does not know how to handle cross sectional dataset, but it "
            "can convert it to a Clinica compliant form (using session ses-M00)"
        )
        super(ClinicaException, self).__init__(self.msg)
