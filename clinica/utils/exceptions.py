"""This module handles Clinica exceptions."""


class ClinicaException(Exception):
    """Base class for Clinica exceptions."""


class ClinicaMissingDependencyError(ClinicaException):
    """Base class for Clinica dependencies errors."""


class ClinicaBIDSError(ClinicaException):
    """Base class for BIDS errors."""


class ClinicaCAPSError(ClinicaException):
    """Base class for CAPS errors."""


class ClinicaParserError(ClinicaException):
    """Base class for parser errors."""


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
