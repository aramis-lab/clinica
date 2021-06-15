__all__ = ["__version__"]

import pkgutil

# Load the Clinica package version
import sys

__version__ = pkgutil.get_data(__package__, "VERSION").decode("ascii").strip()
version = __version__

# import pkg_resources
# version = pkg_resources.require("Clinica")[0].version
# __version__ = version

# python 3.6 minimum version is required
if sys.version_info < (3, 6):
    print(f"Clinica {__version__} requires Python 3.6")
    sys.exit(1)

# Note: The following lines have been commented for speed purposes.
# from clinica.utils.stream import cprint
# import __builtin__
# __builtin__.cprint = cprint
