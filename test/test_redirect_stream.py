from nose.tools import assert_equals
from clinica.utils.redirect_stream import RedirectStream

def test_RedirectStreamPrint():
    stream = RedirectStream()
    stream.active()
    print "hey"
    stream.desactive()
    assert_equals(stream.get_stdout(), "hey\n")
