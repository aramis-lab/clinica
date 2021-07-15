from click_option_group import optgroup

option = optgroup.option

pipeline_options = optgroup.group("Pipeline options")

standard_options = optgroup.group("Standard options")

advanced_options = optgroup.group("Advanced options")


def custom_pipeline_options(name: str, *args, **kwargs):
    return optgroup.group(name=name, *args, **kwargs)
