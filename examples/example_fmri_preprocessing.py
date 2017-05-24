#!/usr/bin/python#
# -*- coding: utf-8 -*-

from pipelines.fmri_preprocessing.fmri_preprocessing_pipeline import fMRIPreprocessing

pipeline = fMRIPreprocessing('/Users/jeremy.guillon/Repositories/multiconproject/data/HMTC_BIDS',
                             '/Users/jeremy.guillon/Repositories/multiconproject/data/HMTC_CAPS')
pipeline.parameters = {
    'num_slices': 45,
    'time_repetition': 2.4,
    'echo_times': [5.19, 7.65],
    'blip_direction': 1,
    'total_readout_time': 15.6799
}

pipeline.base_dir = '/Users/jeremy.guillon/Tmp2'
pipeline.build_core_nodes()
pipeline.build_input_node()
pipeline.input_node._check_inputs(pipeline.get_input_fields())
pipeline.input_node._got_inputs
# pipeline.write_graph()

pipeline.run(plugin='MultiProc', plugin_args={'n_procs' : 2})
