# coding: utf8

"""
   Test to verify coding style dans clinica
"""

__author__ = "Mauricio Diaz"
__copyright__ = "Copyright 2016-2019 The Aramis Lab Team"
__credits__ = ["Mauricio Diaz"]
__license__ = "See LICENSE.txt file"
__version__ = "0.2.0"
__maintainer__ = "Mauricio Diaz"
__email__ = "mauricio.diaz@inria.fr"
__status__ = "Development"

import pycodestyle


def test_coding_style():
    """Test that we conform to PEP-8."""
    style = pycodestyle.StyleGuide(
            quiet=False,
            ignore=['E203', 'E121', 'E123', 'E126', 'E133',
                    'E226', 'E241', 'E242', 'E704', 'W503',
                    'E501', 'W504', 'W505', 'W605'])
    result = style.check_files(['clinica/'])
    result.print_statistics()
    assert result.total_errors ==  0, "Found code style errors (and warnings)."
    pass

