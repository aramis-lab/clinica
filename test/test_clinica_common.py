from nose.tools import assert_equals

from clinica.cmdline import PipelineLoader

def test_PiepelineLoader():
    import os
    env_good = os.getcwd() + "/data/emulate_home";
    env = env_good+ ":/not-existing/dir:/another-not-existing-dir:..w.r.o.n.g...stuff: -- ::"
    pipeline_src1 = env_good+ "/pipeline_1/src"
    pipeline_src2 = env_good+ "/pipeline_2/src"
    os.environ['TEST_CLINICAPATH'] = env
    loader = PipelineLoader(env)

    assert_equals(loader.extract_existing_paths(env), [env_good])
    assert_equals(loader.discover_path_with_subdir([env_good], 'src'), [pipeline_src1, pipeline_src2])