# coding: utf8

from os.path import dirname
from os.path import abspath
from os.path import join as pathjoin

from grabbit import Layout

__all__ = ['CAPSLayout']


class CAPSLayout(Layout):

    def __init__(self, path, config=None, **kwargs):
        if config is None:
            root = dirname(abspath(__file__))
            config = pathjoin(root, 'config', 'caps.json')

        super(CAPSLayout, self).__init__(path, config,
                                         dynamic_getters=True, **kwargs)
