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
