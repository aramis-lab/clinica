from nose.tools import assert_equals, assert_true

from clinica.cmdline import PipelineLoader

def test_PiepelineLoader():
    from clinica.pipeline import Pipeline
    from test.data.emulate_home.pipeline_1.src.pipeline_1_cli import TestPip1
    from test.data.emulate_home.base_class import base
    import os
    env_good = os.getcwd() + "/data/emulate_home";
    env = env_good+ ":/not-existing/dir:/another-not-existing-dir:..w.r.o.n.g...stuff: -- ::"
    pipeline_src1 = env_good+ "/pipeline_1/src"
    pipeline_file1 = env_good+ "/pipeline_1/src/pipeline_1_cli.py"
    pipeline_src2 = env_good+ "/pipeline_2/src"
    pipeline_file2 = env_good + "/pipeline_2/src/pipeline_2_cli.py"
    env_name = 'TEST_CLINICAPATH'
    os.environ[env_name] = env
    loader = PipelineLoader(env=env_name,baseclass=base)

    assert_equals(loader.extract_existing_paths(env), [env_good])
    assert_equals(loader.discover_path_with_subdir([env_good], 'src'), [pipeline_src1, pipeline_src2])
    assert_equals(loader.find_files([pipeline_src1, pipeline_src2], r".*_cli\.py$"), [pipeline_file1, pipeline_file2])
    assert_equals(loader.load_class(base, pipeline_file1).ret(), 1)
    assert_equals(loader.load_class(base, pipeline_file2).ret(), 2)
    classes = loader.load()
    assert_equals(len(classes), 2)
    assert_equals(classes[0].ret(), 1)
    assert_equals(classes[1].ret(), 2)

