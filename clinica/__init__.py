import os
from importlib.metadata import version

__all__ = ["__version__"]

__version__ = version("clinica")

# Disable Nipype / Pydra telemetry.
os.environ["NO_ET"] = "yes"
