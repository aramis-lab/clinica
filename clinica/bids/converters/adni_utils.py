
def replace_sequence_chars(sequence_name):
    import re
    return re.sub('[ /;*():]', '_', sequence_name)


def fill_zeros(s, length):
    return ('0' * (length - len(str(s)))) + str(s)


def days_between(d1, d2):
    '''
    Calculate the days between two dates

    :param d1: date 1
    :param d2: date 2
    :return: number of days between date 2 and date 1
    '''
    from datetime import datetime
    d1 = datetime.strptime(d1, "%Y-%m-%d")
    d2 = datetime.strptime(d2, "%Y-%m-%d")
    return abs((d2 - d1).days)


def viscode_to_session(viscode):
    if viscode == 'bl':
        return 'M00'
    else:
        return viscode.capitalize()


def center_nifti_origin(input_image, output_image):
    '''
     Put the origin of the coordinate system at the center of the image.

    :param input_image: path to the input image
    :param output_image: path to the output image (where the result will be stored)
    :return:
    '''

    import nibabel as nib
    import numpy as np

    img = nib.load(input_image)
    canonical_img = nib.as_closest_canonical(img)
    hd = canonical_img.header
    if hd['quatern_b'] != 0 or hd['quatern_c'] != 0 or hd['quatern_d'] != 0:
        print 'Warning: Not all values in quatern are equal to zero'
    qform = np.zeros((4, 4))
    for i in range(1, 4):
        qform[i - 1, i - 1] = hd['pixdim'][i]
        qform[i - 1, 3] = -1.0 * hd['pixdim'][i] * hd['dim'][i] / 2.0
    new_img = nib.Nifti1Image(canonical_img.get_data(caching='unchanged'), qform)
    nib.save(new_img, output_image)






