from nose.tools import assert_equals

from clinica.cmdline import PipelineLoader

def test_PiepelineLoader():
    import os
    env_good = os.getcwd() + "/data/emulate_home";
    env = env_good+ ":/not-existing/dir:/another-not-existing-dir:..w.r.o.n.g...stuff: -- ::"
    os.environ['TEST_CLINICAPATH'] = env
    loader = PipelineLoader(env)

    assert_equals(loader.extract_paths(env), env)