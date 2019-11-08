# coding: utf8


def get_pipeline_parameters(full_width_at_half_maximum=None,
                            t1_native_space=None,
                            freesurfer_brain_mask=None,
                            unwarping=None):
    parameters = {
        'full_width_at_half_maximum': full_width_at_half_maximum or [8, 8, 8],
        't1_native_space': t1_native_space or False,
        'freesurfer_brain_mask': freesurfer_brain_mask or False,
        'unwarping': unwarping or False,
    }

    return parameters
