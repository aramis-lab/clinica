def prepare_flowfields(flow_fields, tissues):
    return [[f] * len(tissues) for f in flow_fields]


def join_smoothed_files(smoothed_normalized_files):
    """Join outputs."""
    return [
        [x for smooth in subject for x in smooth]
        for subject in zip(*smoothed_normalized_files)
    ]
