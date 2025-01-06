import os
from importlib.metadata import version

from packaging.version import Version

__all__ = ["__version__", "get_version"]


__version__ = version("clinica")


def get_version() -> Version:
    return Version(__version__)


# Disable Nipype / Pydra telemetry.
os.environ["NO_ET"] = "yes"
