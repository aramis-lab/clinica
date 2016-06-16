from __future__ import print_function
import sys
from optparse import OptionParser

import clinica



def execute():
    """
    clinica <command> [path=current_directory] [options]
    ex:
    $cd WorkingDir/
    $clinica visualize --id=1,2,3
    """
    parser = OptionParser()
    parser.add_option("-i", "--id", dest="id",
                      help="unique identifier")
    parser.add_option("-q", "--quiet",
                      action="store_false", dest="verbose", default=True,
                      help="don't print status messages to stdout")

    (options, args) = parser.parse_args()
    print(args)
    print(options)

if __name__ == '__main__':
    execute()
