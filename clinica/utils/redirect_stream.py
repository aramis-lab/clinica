"""
Redirect stream and log
"""
import sys
from cStringIO import StringIO

class RedirectStream:

    def active(self):
        sys.stdout = self.stdout = StringIO()
        sys.stderr = self.stderr = StringIO()

    def desactive(self):
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__

    def get_stdout(self):
        return self.stdout.getvalue()

    def get_stderr(self):
        return self.stderr.getvalue()


def Message(msg):pass
