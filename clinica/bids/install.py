from __future__ import print_function
from os.path import basename, realpath,split,join
from shutil import copyfile

def add_bids_converter(src):
    copy_data(src, 'converters')

def add_bids_data(src):
    copy_data(src, 'data')

def copy_data(src, dest):
    converter_path = join(split(realpath(__file__))[0], dest)
    src_name = basename(src)
    dst = join(converter_path, src_name)
    try:
        copyfile(src, dst)
    except:
        print('re-run with sudo permision')
        exit(-1)

    print('Installation of %s done' % basename(src))

