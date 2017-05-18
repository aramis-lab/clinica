#!/usr/bin/python#
# -*- coding: utf-8 -*-

"""

"""

from clinica.pipeline.fmri.preprocessing import FMRIPreprocessing
import nipype.pipeline.engine as npe
import clinica.pipeline.engine as cpe


# Pipeline 1
# ==========
pipeline1 = FMRIPreprocessing(name="A", bids_dir='/Users/jeremy.guillon/Repositories/multiconproject/data/HMTC_BIDS')
pipeline1.parameters = {
    'num_slices': 45,
    'time_repetition': 2.4,
    'echo_times': [5.19, 7.65],
    'blip_direction': 1,
    'total_readout_time': 15.6799
}
pipeline1.build_core_nodes()
pipeline1.build_input_node()


# Pipeline 2
# ==========
pipeline2 = FMRIPreprocessing(name="B", caps_dir='/Users/jeremy.guillon/Repositories/multiconproject/data/HMTC_CAPS')
pipeline2.parameters = {
    'num_slices': 45,
    'time_repetition': 2.4,
    'echo_times': [5.19, 7.65],
    'blip_direction': 1,
    'total_readout_time': 15.6799
}
pipeline2.build_core_nodes()
pipeline2.build_output_node()


# Workflow of pipelines
# =====================
wf = npe.Workflow("BigOne")
wf.connect(pipeline1, 'Output.mc_params', pipeline2, 'Input.phasediff')
wf.base_dir = '/Users/jeremy.guillon/Tmp2'
wf.write_graph()