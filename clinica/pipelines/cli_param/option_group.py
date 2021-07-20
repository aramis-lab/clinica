from click_option_group import optgroup

option = optgroup.option

pipeline_specific_options = optgroup.group(
    "Pipeline-specific options", help="Options specific to the pipeline being run"
)

common_pipelines_options = optgroup.group(
    "Common pipelines options", help="Options common to all Clinica pipelines"
)

advanced_pipeline_options = optgroup.group(
    "Advanced pipeline options", help="For experts only"
)


def custom_pipeline_options(name: str, *args, **kwargs):
    return optgroup.group(name=name, *args, **kwargs)
