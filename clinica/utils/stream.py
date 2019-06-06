# coding: utf8

"""
Redirect stream and log
"""
import sys

clinica_verbose = False


class FilterOut(object):
    def __init__(self, stdout):
        self.stdout = stdout

    def write(self, text):
        import re
        if not text:
            return
        if re.match('^(@clinica@)', text):
            self.stdout.write(text.replace("@clinica@", ""))
            self.stdout.flush()

    def flush(self): self.stdout.flush()

    def __enter__(self):
        self.origin_stdout = sys.stdout
        self.flush()

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout = self.origin_stdout
        self.flush()


def active_cprint():
    sys.stdout = FilterOut(sys.stdout)


def cprint(msg):
    global clinica_verbose
    if clinica_verbose is True:
        print(msg)
    else:
        print("@clinica@ %s" % msg)
    sys.stdout.flush()
