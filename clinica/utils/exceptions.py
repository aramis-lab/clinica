"""

"""

class ClinicaException(Exception):
    """ Base class for Clinica exceptions """
    def __init__(self, *args, **kwargs):
        pass

class ClinicaMissingDependencyError(ClinicaException):
    """ Base class for Clinica dependencies errors """
    def __init__(self, *args, **kwargs):
        pass

class ClinicaBIDSError(ClinicaException):
    """ Base class for BIDS errors """
    def __init__(self, *args, **kwargs):
        pass

class ClinicaCAPSError(ClinicaException):
    """ Base class for CAPS errors """
    def __init__(self, *args, **kwargs):
        pass