import numpy as np
import pandas as pd
import nibabel as nib
import csv
from ntpath import basename, splitext
from scipy.spatial.distance import squareform
from os.path import join
import scipy.io
from skimage import img_as_float


def get_caps_t1_list(input_directory, subjects_visits_tsv, atlas_id):
    ''''
    path to arrive to the list of the file with the statistics on atlas_id 
    '''

    subjects_visits = pd.io.parsers.read_csv(subjects_visits_tsv, sep='\t')
    if list(subjects_visits.columns.values) != ['participant_id', 'session_id']:
        raise Exception('Subjects and visits file is not in the correct format.')
    subjects = list(subjects_visits.participant_id)
    sessions = list(subjects_visits.session_id)
    image_list = [join(input_directory, '/subjects/' + subjects[i] + '/'
                       + sessions[i] + '/t1/spm/atlas_statistics/'+subjects[i]+'_'+sessions[i]+'_T1wspace-'+atlas_id+'_map-graymatter_statistics.tsv') for i in range(len(subjects))]
    return image_list

def get_caps_pet_list(input_directory, subjects_visits_tsv, group_id, atlas_id):

    subjects_visits = pd.io.parsers.read_csv(subjects_visits_tsv, sep='\t')
    if list(subjects_visits.columns.values) != ['participant_id', 'session_id']:
        raise Exception('Subjects and visits file is not in the correct format.')
    subjects = list(subjects_visits.participant_id)
    sessions = list(subjects_visits.session_id)
    image_list = [join(input_directory, 'analysis-series-default/subjects/' + subjects[i] + '/'
                       + sessions[i] + '/pet/atlas_statistics/'+ subjects[i]+'_'+sessions[i]+'_space-'+atlas_id+'_map-fdgstatistic2.tsv') for i in range(len(subjects))]
    return image_list

def load_data(image_list,subjects_visits_tsv):
    subjects_visits = pd.io.parsers.read_csv(subjects_visits_tsv, sep='\t')
    subjects = list(subjects_visits.participant_id)
    subj_average=[]
    all_vector=np.array([])
    data_=np.array([])
    read_file= pd.io.parsers.read_csv(image_list[0], sep='\t',usecols=[2],header=0)
    read_file=read_file.mean_scalar
    data=np.zeros((len(subjects), len(read_file)))
    for i in xrange(len(image_list)):
        tsv_file= pd.io.parsers.read_csv(image_list[i], sep='\t',usecols=[2],header=0)
        subj_average=tsv_file.mean_scalar
        all_vector=np.append(all_vector,subj_average)
    data_=np.split(all_vector, len(image_list))
    for i in xrange(len(image_list)):
        for j in xrange(len(subj_average)):
            data[i][j]=data_[i][j]
    scipy.io.savemat('data.mat', {'data':data})
    return data

def features_weights(image_list, dual_coefficients, sv_indices, scaler=None, mask=None):

    if len(sv_indices) != len(dual_coefficients):
        print "Length dual coefficients: " + str(len(dual_coefficients))
        print "Length indices: " + str(len(sv_indices))
        raise ValueError('The number of support vectors indices and the number of coefficients must be the same.')

    if len(image_list) == 0:
        raise ValueError('The number of images must be greater than 0.')

    sv_images = [image_list[i] for i in sv_indices]

    shape = pd.io.parsers.read_csv(sv_images[0], sep='\t',usecols=[2],header=0)
    weights = np.zeros(len(shape))

    for i in range(len(sv_images)):
        subj = pd.io.parsers.read_csv(sv_images[i], sep='\t',usecols=[2],header=0)
        subj_data = subj.mean_scalar
        weights += dual_coefficients[i] * subj_data

    return weights

def weights_to_nifti(input_image_atlas,weights):
    #input_image_atlas = os.path.join('....', 'atlas_id.nii')

    atlas_image = nib.load(input_image_atlas).get_data()
    labels = list(set(atlas_image.ravel()))
    output_image_weights=np.array(atlas_image,dtype='f')
    for i, n in enumerate(labels):
        index = np.array(np.where(atlas_image== n))
        output_image_weights[index[0, :], index[1, :], index[2, :]]=weights[i]
    output_image = nib.Nifti1Image(output_image_weights, nib.load(input_image_atlas).get_affine(),nib.load(input_image_atlas).get_header())

    return output_image
