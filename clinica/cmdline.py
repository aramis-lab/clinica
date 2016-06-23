from __future__ import print_function
from optparse import OptionParser
from clinica.engine.cworkflow import *
import sys
import os
import subprocess

def visualize(clinicaWorkflow, ids, rebase=False):
    if not clinicaWorkflow.data.has_key('visualize'):
        print("No visualization was defined")
        exit(0)

    class chdir:
        def __init__(self, base):
            self.pwd = os.getcwd()
            os.chdir(base)
        def __del__(self):
            os.chdir(self.pwd)

    change_directory = None
    if rebase is False:
        change_directory = chdir(clinicaWorkflow.base_dir)
    else:
        change_directory = chdir(rebase)

    program, arguments, matchs = clinicaWorkflow.data['visualize']

    def run_program(id): subprocess.Popen([program] + arguments.replace("${%s}" % matchs, id).strip().split(" "))
    [run_program(id) for id in ids]

def shell(clinicaWorkflow):
    workflow = clinicaWorkflow
    __banner__ = "workflow variable is instantiated for you!"
    namespace = globals().copy()
    namespace.update(locals())

    def load_python_shell():
        import readline
        import code
        shell = code.InteractiveConsole(namespace)
        shell.interact(banner=__banner__)

    def load_ipython_shell():
        from IPython.terminal.embed import InteractiveShellEmbed
        InteractiveShellEmbed(user_ns=namespace,banner1=__banner__)()

    try:
        load_ipython_shell()
    except:
        try:
            load_python_shell()
        except:
            print("Impossible to load ipython or python shell")

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
    parser.add_option("-r", "--rebase", dest="rebase",
                      default=False,
                      help="unique identifier")
    parser.add_option("-q", "--quiet",
                      action="store_false", dest="verbose", default=True,
                      help="don't print status messages to stdout")

    (options, args) = parser.parse_args()

    if args[0] == 'visualize':
        if options.id is False:
            print("Missing --id")
            exit(0)
        visualize(load_conf(args[1:]), options.id.split(","), options.rebase)

    if args[0] == 'shell':
        shell(load_conf(args[1:]))


if __name__ == '__main__':
    execute()
