# coding: utf8


def check_group_label(group_label):
    """Check that `group_label` is compliant with specifications."""
    if not group_label.isalnum():
        raise ValueError('Not valid group_label value: it must be composed only by letters and/or numbers '
                         '(given value: %s).' % group_label)


def extract_group_ids(caps_directory):
    """Extract list of group IDs (e.g. ['group-AD', 'group-HC']) based on `caps_directory`/groups folder."""
    import os

    try:
        group_ids = os.listdir(os.path.join(caps_directory, 'groups'))
    except FileNotFoundError:
        group_ids = ['']

    return group_ids
