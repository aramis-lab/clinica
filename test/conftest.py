
"""
    This file contains a set of functional tests designed to check the correct execution of the pipeline and the
    different functions available in Clinica
"""

__author__ = "Mauricio Diaz"
__copyright__ = "Copyright 2016-2019 The Aramis Lab Team"
__credits__ = ["Mauricio Diaz"]
__license__ = "See LICENSE.txt file"
__version__ = "0.2.0"
__maintainer__ = "Mauricio Diaz"
__email__ = "mauricio.diaz@inria.fr"
__status__ = "Development"


import pytest

def pytest_addoption(parser):
    parser.addoption(
            "--working_directory", 
        action="store", 
        help="Working directory for tests"
    )


@pytest.fixture
def cmdopt(request):
    return request.config.getoption("--working_directory")

