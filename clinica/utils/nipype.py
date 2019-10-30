# coding: utf8


def fix_join(path, *paths):
    # This workaround is used in pipelines like DWIPreprocessingUsingT1
    # In the workflow.connect part, you can use some function that are used as string, causing an import error
    import os
    return os.path.join(path, *paths)
