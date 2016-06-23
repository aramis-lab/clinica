from __future__ import print_function
import sys
from optparse import OptionParser
from clinica.engine.cworkflow import *
import os
import subprocess

def visualize(clinicaWorkflow, ids):
    if not clinicaWorkflow.data.has_key('visualize'):
        print("No visualization was defined")
        exit(0)

    class chdir:
        def __init__(self, base):
            self.pwd = os.getcwd()
            os.chdir(base)
        def __del__(self):
            os.chdir(self.pwd)

    change_directory = chdir(clinicaWorkflow.base_dir)
    program, arguments, matchs = clinicaWorkflow.data['visualize']

    def run_program(id): subprocess.Popen([program] + arguments.replace("${%s}" % matchs, id).strip().split(" "))
    [run_program(id) for id in ids]

    print("exec %s %s %s" % (program, arguments, matchs))

def load_conf(args):
    import cPickle

    def load(path):
        file = os.path.join(path, "clinica.pkl")
        if os.path.isfile(file): return cPickle.load(open(file))
        return False

    wk = False

    if len(args) == 0:
        wk = load(os.getcwd())
    elif os.path.isdir(args[0]):
        wk = load(args[0])

    if not wk:
        print("No <clinica.pkl> file found!")
        exit(0)

    return wk

def execute():
    """
    clinica <command> [path=current_directory] [options]
    ex:
    $cd WorkingDir/
    $clinica visualize --id=1,2,3
    """
    parser = OptionParser()
    parser.add_option("-i", "--id", dest="id",
                      default=False,
                      help="unique identifier")
    parser.add_option("-q", "--quiet",
                      action="store_false", dest="verbose", default=True,
                      help="don't print status messages to stdout")

    (options, args) = parser.parse_args()
    print(args)
    print(options)
    if args[0] == 'visualize':
        if options.id is False:
            print("Missing --id")
            exit(0)
        visualize(load_conf(args[1:]), options.id.split(","))
        print('visualize....')

if __name__ == '__main__':
    execute()
