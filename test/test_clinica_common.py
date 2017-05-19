from nose.tools import assert_equals

from clinica.cmdline import ClinicaClassLoader

def test_PiepelineLoader():
    from test.data.emulate_home.base_class import base
    import os
    env_good = os.getcwd() + "/data/emulate_home";
    env_pipeline = env_good + "/pipelines"
    pipeline_src1 = env_pipeline+ "/pipeline_1"
    pipeline_file1 = env_pipeline+ "/pipeline_1/pipeline_1_cli.py"
    pipeline_src2 = env_pipeline+ "/pipeline_2"
    pipeline_file2 = env_pipeline + "/pipeline_2/pipeline_2_cli.py"
    env_name = 'TEST_CLINICAPATH'
    os.environ[env_name] = env_good
    loader = ClinicaClassLoader(env=env_name, baseclass=base, extra_dir="pipelines")

    assert_equals(loader.discover_path_with_subdir(env_pipeline), [pipeline_src1, pipeline_src2])
    assert_equals(loader.find_files([pipeline_src1, pipeline_src2], r".*_cli\.py$"), [pipeline_file1, pipeline_file2])
    assert_equals(loader.load_class(base, pipeline_file1).ret(), 1)
    assert_equals(loader.load_class(base, pipeline_file2).ret(), 2)
    classes = loader.load()
    assert_equals(len(classes), 2)
    assert_equals(classes[0].ret(), 1)
    assert_equals(classes[1].ret(), 2)

def test_ConvertersLoader():
    from test.data.emulate_home.base_class import base
    import os
    env_good = os.getcwd() + "/data/emulate_home";
    env_converter = env_good + "/iotools/converters"
    env_name = 'TEST_CLINICAPATH'
    os.environ[env_name] = env_good
    loader = ClinicaClassLoader(env=env_name, baseclass=base, extra_dir="iotools/converters", reg=r".*_bids\.py$")
    classes = loader.load()
    assert_equals(len(classes), 1)
    assert_equals(classes[0].ret(), "x_to_bids")

