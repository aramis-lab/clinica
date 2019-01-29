# coding: utf8


import numpy as np
import pandas as pd
import nibabel as nib
from os.path import join

__author__ = "Jorge Samper-Gonzalez"
__copyright__ = "Copyright 2016-2019 The Aramis Lab Team"
__credits__ = ["Jorge Samper-Gonzalez", "Simona Bottani"]
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__maintainer__ = "Jorge Samper-Gonzalez"
__email__ = "jorge.samper-gonzalez@inria.fr"
__status__ = "Development"


def get_caps_t1_list(input_directory, subjects_visits_tsv, group_id, atlas_id):
    """
    path to arrive to the list of the file with the statistics on atlas_id
    Args:
        input_directory:
        subjects_visits_tsv:
        group_id:
        atlas_id:

    Returns:

    """

    from os.path import join
    import pandas as pd

    subjects_visits = pd.io.parsers.read_csv(subjects_visits_tsv, sep='\t')
    if list(subjects_visits.columns.values) != ['participant_id', 'session_id']:
        raise Exception('Subjects and visits file is not in the correct format.')
    subjects = list(subjects_visits.participant_id)
    sessions = list(subjects_visits.session_id)
    image_list = [join(input_directory + '/subjects/' + subjects[i] + '/'
                       + sessions[i] + '/t1/spm/dartel/group-' + group_id + '/atlas_statistics/' + subjects[i] + '_'
                       + sessions[i]+'_T1w_space-'+atlas_id+'_map-graymatter_statistics.tsv')
                  for i in range(len(subjects))]
    return image_list


def get_caps_pet_list(input_directory, subjects_visits_tsv, group_id, atlas_id):
    """

    Args:
        input_directory:
        subjects_visits_tsv:
        group_id:
        atlas_id:

    Returns:

    """

    subjects_visits = pd.io.parsers.read_csv(subjects_visits_tsv, sep='\t')
    if list(subjects_visits.columns.values) != ['participant_id', 'session_id']:
        raise Exception('Subjects and visits file is not in the correct format.')
    subjects = list(subjects_visits.participant_id)
    sessions = list(subjects_visits.session_id)
    image_list = [join(input_directory, 'analysis-series-default/subjects/' + subjects[i] + '/'
                       + sessions[i] + '/pet/atlas_statistics/' + subjects[i] + '_' + sessions[i]
                       + '_space-' + atlas_id + '_map-fdgstatistic2.tsv')
                  for i in range(len(subjects))]
    return image_list


def load_data(image_list, subjects):
    """

    Args:
        image_list:
        subjects:

    Returns:

    """

    subj_average = []
    all_vector = np.array([])
    read_file = pd.io.parsers.read_csv(image_list[0], sep='\t', usecols=[2], header=0)
    read_file = read_file.mean_scalar
    data = np.zeros((len(subjects), len(read_file)))

    for i in range(len(image_list)):
        tsv_file = pd.io.parsers.read_csv(image_list[i], sep='\t', usecols=[2], header=0)
        subj_average = tsv_file.mean_scalar
        all_vector = np.append(all_vector, subj_average)
    data_temp = np.split(all_vector, len(image_list))

    for i in range(len(image_list)):
        for j in range(len(subj_average)):
            data[i][j] = data_temp[i][j]
    return data


def features_weights(image_list, dual_coefficients, sv_indices, scaler=None):
    """

    Args:
        image_list:
        dual_coefficients:
        sv_indices:
        scaler:

    Returns:

    """

    if len(sv_indices) != len(dual_coefficients):
        print("Length dual coefficients: " + str(len(dual_coefficients)))
        print("Length indices: " + str(len(sv_indices)))
        raise ValueError('The number of support vectors indices and the number of coefficients must be the same.')

    if len(image_list) == 0:
        raise ValueError('The number of images must be greater than 0.')

    sv_images = [image_list[i] for i in sv_indices]

    shape = pd.io.parsers.read_csv(sv_images[0], sep='\t', usecols=[2], header=0)
    weights = np.zeros(len(shape))

    for i in range(len(sv_images)):
        subj = pd.io.parsers.read_csv(sv_images[i], sep='\t', usecols=[2], header=0)
        subj_data = subj.mean_scalar
        weights += dual_coefficients[i] * subj_data

    return weights


def weights_to_nifti(weights, atlas, output_filename):
    """

    Args:
        atlas:
        weights:
        output_filename:

    Returns:

    """

    from os.path import join, split, realpath

    from clinica.utils.atlas import AtlasAbstract

    atlas_path = None
    atlas_classes = AtlasAbstract.__subclasses__()
    for atlas_class in atlas_classes:
        if atlas_class.get_name_atlas() == atlas:
            atlas_path = atlas_class.get_atlas_labels()

    if not atlas_path:
        raise ValueError('Atlas path not found for atlas name ' + atlas)

    atlas_image = nib.load(atlas_path)
    atlas_data = atlas_image.get_data()
    labels = list(set(atlas_data.ravel()))
    output_image_weights = np.array(atlas_data, dtype='f')
    for i, n in enumerate(labels):
        index = np.array(np.where(atlas_data == n))
        output_image_weights[index[0, :], index[1, :], index[2, :]] = weights[i]

    affine = atlas_image.get_affine()
    header = atlas_image.get_header()
    output_image = nib.Nifti1Image(output_image_weights, affine, header)
    nib.save(output_image, output_filename)
