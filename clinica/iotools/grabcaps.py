# coding: utf8

from os.path import dirname
from os.path import realpath
from os.path import join as pathjoin

from grabbit import Layout

__all__ = ['CAPSLayout']


class CAPSLayout(Layout):

    def __init__(self, path, config=None, **kwargs):
        import os
        if config is None:
            root = dirname(realpath(__file__))
            config = pathjoin(os.path.abspath(root), '..', 'resources', 'layouts', 'caps.json')
        super(CAPSLayout, self).__init__(path, config,
                                         dynamic_getters=True, **kwargs)


# if __name__ == '__main__':
#     layout = CAPSLayout(path='/Users/jeremy.guillon/Repositories/multiconproject/data/HMTC_CAPS')
#     #print(layout.get(target='subject', return_type='dir', regex_search=True))
#     print(layout.get(target='session', return_type='dir', regex_search=False))
#     print(layout.get(fmripreprocessing='bold_motionparams'))
#     print(layout.get(target='modality', return_type='dir', regex_search=False))
#     print(layout.get(target='modality', return_type='id', regex_search=False))
