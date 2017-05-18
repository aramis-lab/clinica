#!/usr/bin/python#
# -*- coding: utf-8 -*-

from clinica.pipeline.fmri.preprocessing import FMRIPreprocessing

pipeline = FMRIPreprocessing('/Users/jeremy.guillon/Repositories/multiconproject/data/HMTC_BIDS',
                             '/Users/jeremy.guillon/Repositories/multiconproject/data/HMTC_CAPS')
pipeline.parameters = {
    'num_slices': 45,
    'time_repetition': 2.4,
    'echo_times': [5.19, 7.65],
    'blip_direction': 1,
    'total_readout_time': 15.6799
}

pipeline.base_dir = '/Users/jeremy.guillon/Tmp2'
pipeline.build()
pipeline.write_graph()

pipeline.run(plugin='MultiProc', plugin_args={'n_procs' : 2})
