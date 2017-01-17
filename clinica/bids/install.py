from __future__ import print_function
from os.path import basename, realpath,split,join
from shutil import copyfile

def add_bids_converter(src):
    converter_path = join(split(realpath(__file__))[0], 'converters')
    src_name = basename(src)
    dst = join(converter_path, src_name)
    try:
        copyfile(src, dst)
    except PermissionError:
        print('re-run with sudo permision')

    print('copy done')
