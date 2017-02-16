#!/usr/bin/python#
# -*- coding: utf-8 -*-

from clinica.pipeline.fmri.preprocessing import FMRIPreprocessing

pipeline = FMRIPreprocessing('/Volumes/dataARAMIS/users/jeremy.guillon/multiconproject/data/HMTC_BIDS',
                             '/Volumes/dataARAMIS/users/jeremy.guillon/multiconproject/data/HMTC_CAPS')
pipeline.parameters = {
    'num_slices': 45,
    'time_repetition': 2.4
}
pipeline.base_dir = '/Users/jeremy.guillon/Tmp'
pipeline.run()
