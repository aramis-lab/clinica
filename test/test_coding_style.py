"""
   Test tto verify coding style dans clinica
"""

__author__ = "Mauricio Diaz"
__copyright__ = "Copyright 2016-2019 The Aramis Lab Team"
__credits__ = ["Arnaud Marcoux"]
__license__ = "See LICENSE.txt file"
__version__ = "0.2.0"
__maintainer__ = "Mauricio Diaz"
__email__ = "mauricio.diaz@inria.fr"
__status__ = "Development"


import warnings
import sys
import pycodestyle


def test_coding_style():
    """Test that we conform to PEP-8."""
    style = pycodestyle.StyleGuide(quiet=True)
    result = style.check_files(['file1.py', 'file2.py'])
    self.assertEqual(result.total_errors, 0,
                     "Found code style errors (and warnings).")    
    pass

