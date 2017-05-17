#!/usr/bin/python#
# -*- coding: utf-8 -*-

"""

"""

from clinica.pipeline.fmri.preprocessing import FMRIPreprocessing

pipeline1 = FMRIPreprocessing('/Users/jeremy.guillon/Repositories/multiconproject/data/HMTC_BIDS',
                              '/Users/jeremy.guillon/Repositories/multiconproject/data/HMTC_CAPS', name="A")
pipeline1.parameters = {
    'num_slices': 45,
    'time_repetition': 2.4,
    'echo_times': [5.19, 7.65],
    'blip_direction': 1,
    'total_readout_time': 15.6799
}
pipeline1.base_dir = '/Users/jeremy.guillon/Tmp2'
pipeline1.build_core_nodes()
pipeline1.build_input_node()

pipeline2 = FMRIPreprocessing('/Users/jeremy.guillon/Repositories/multiconproject/data/HMTC_BIDS',
                              '/Users/jeremy.guillon/Repositories/multiconproject/data/HMTC_CAPS', name="B")
pipeline2.parameters = {
    'num_slices': 45,
    'time_repetition': 2.4,
    'echo_times': [5.19, 7.65],
    'blip_direction': 1,
    'total_readout_time': 15.6799
}
pipeline2.base_dir = '/Users/jeremy.guillon/Tmp2'
pipeline2.build_core_nodes()
pipeline2.build_output_node()

# clinica.engine.connect([pipeline1:'Output.mc_params', pipeline2:'Input.phasediff'])



import nipype.pipeline.engine as npe

wf = npe.Workflow("BigOne")
wf.connect(pipeline1, 'Output.mc_params', pipeline2, 'Input.phasediff')
wf.write_graph()