

def zip_nii(in_file):
    from os import getcwd
    from os.path import abspath, join
    import gzip
    import shutil
    from nipype.utils.filemanip import split_filename

    _, base, ext = split_filename(in_file)

    # Already compressed
    if ext[-3:].lower() == ".gz":
        return in_file
    # Not compressed
    out_file = abspath(join(getcwd(), base + ext + '.gz'))
    with open(in_file, 'rb') as f_in, gzip.open(out_file, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)

    return out_file


def unzip_nii(in_file):
    from nipype.utils.filemanip import split_filename
    from nipype.algorithms.misc import Gunzip

    _, base, ext = split_filename(in_file)

    # Not compressed
    if ext[-3:].lower() != ".gz":
        return in_file
    # Compressed
    gunzip = Gunzip(in_file=in_file)
    gunzip.run()
    return gunzip.aggregate_outputs().out_file


