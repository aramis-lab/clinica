import numpy as np
import pandas as pd
import nibabel as nib
import csv
from ntpath import basename, splitext
from scipy.spatial.distance import squareform
from os.path import join
import scipy.io
from skimage import img_as_float


def get_caps_t1_list(input_directory, subjects_visits_tsv, group_id, atlas_id):

    subjects_visits = pd.io.parsers.read_csv(subjects_visits_tsv, sep='\t')
    if list(subjects_visits.columns.values) != ['participant_id', 'session_id']:
        raise Exception('Subjects and visits file is not in the correct format.')
    subjects = list(subjects_visits.participant_id)
    sessions = list(subjects_visits.session_id)
    image_list = [join(input_directory, 'analysis-series-default/subjects/' + subjects[i] + '/'
                       + sessions[i] + '/t1/spm/atlas_statistics/'+subjects[i]+'_'+sessions[i]+'_space-'+atlas_id+'_map-gm_statistic.tsv') for i in range(len(subjects))]
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

    input_image_atlas=os.path.join('....','atlas_id.nii')
    atlas_image = nib.load(input_image_atlas).get_data()
    labels = list(set(atlas_image.ravel()))
    output_image_weights=np.array(atlas_image,dtype='f')
    for i, n in enumerate(labels):
        index = np.array(np.where(image == n))
        output_image_weights[index[0, :], index[1, :], index[2, :]]=vettore[i]
    output_image = nib.Nifti1Image(output_image_weights, nib.load(input_image_atlas).get_affine(),nib.load(input_image_atlas).get_header())
    #output_image.to_filename('/Users/simona.bottani/Desktop/Database_60_subjects/SVM_output/aicha/pet/AD_vs_CN_balanced__weights.nii') Posso definirlo dopo
    return output_image

def svm_binary_classification(input_image_atlas,image_list, diagnosis_list, output_directory, kernel_function=None, existing_gram_matrix=None, mask_zeros=True, scale_data=False, balanced=False, outer_folds=10, inner_folds=10, n_threads=10, c_range=np.logspace(-10, 2, 1000), save_gram_matrix=False, save_subject_classification=False, save_dual_coefficients=False, scaler=None, data_mask=None, save_original_weights=False, save_features_image=True):

    if (kernel_function is None and existing_gram_matrix is None) | (kernel_function is not None and existing_gram_matrix is not None):
        raise ValueError('Kernel_function and existing_gram_matrix are mutually exclusive parameters.')

    results = dict()
    dx_filter = np.unique(diagnosis_list)

    if kernel_function is not None:
        print 'Loading ' + str(len(image_list)) + ' subjects'
        x0, orig_shape, data_mask = load_data(image_list, mask=mask_zeros)
        print 'Subjects loaded'
        print 'Calculating Gram matrix'

        if scale_data:
            x_all = scale(x0)
        else:
            x_all = x0

        gram_matrix = kernel_function(x_all)
        print 'Gram matrix calculated'

    if existing_gram_matrix is not None:
        gram_matrix = existing_gram_matrix
        if (gram_matrix.shape[0] != gram_matrix.shape[1]) | (gram_matrix.shape[0] != len(image_list)):
            raise ValueError('The existing Gram matrix must be a square matrix with number of rows and columns equal to the number of images.')

    if save_gram_matrix:
        np.savetxt(join(output_directory, 'gram_matrix.txt'), gram_matrix)

    for i in range(len(dx_filter)):
        for j in range(i + 1, len(dx_filter)):
            dx1 = dx_filter[i]
            dx2 = dx_filter[j]

            ind1 = []
            ind2 = []
            for k in range(len(diagnosis_list)):
                if diagnosis_list[k] == dx1:
                    ind1.append(k)
                if diagnosis_list[k] == dx2:
                    ind2.append(k)

            indices = ind1 + ind2

            current_subjects = [image_list[k] for k in indices]
            current_diagnosis = [diagnosis_list[k] for k in indices]

            y = np.array([0] * len(ind1) + [1] * len(ind2))
            gm = gram_matrix[indices, :][:, indices]

            classification_str = dx1 + '_vs_' + dx2 + ('_balanced' if balanced else '_not_balanced')
            print 'Running ' + dx1 + ' vs ' + dx2 + ' classification'

            y_hat, dual_coefficients, sv_indices, intersect, c = nested_folds(gm, y, c_range, balanced=balanced, outer_folds=outer_folds, inner_folds=inner_folds, n_threads=n_threads)

            evaluation = evaluate_prediction(y, y_hat)

            print '\nTrue positive %0.2f' % len(evaluation['predictions'][0])
            print 'True negative %0.2f' % len(evaluation['predictions'][1])
            print 'False positive %0.2f' % len(evaluation['predictions'][2])
            print 'False negative %0.2f' % len(evaluation['predictions'][3])

            print 'Accuracy %0.2f' % evaluation['accuracy']
            print 'Balanced accuracy %0.2f' % evaluation['balanced_accuracy']
            print 'Sensitivity %0.2f' % evaluation['sensitivity']
            print 'Specificity %0.2f' % evaluation['specificity']
            print 'Positive predictive value %0.2f' % evaluation['ppv']
            print 'Negative predictive value %0.2f \n' % evaluation['npv']

            if save_dual_coefficients:
                np.save(join(output_directory, classification_str + '__dual_coefficients'), dual_coefficients[0])
                np.save(join(output_directory, classification_str + '__sv_indices'), sv_indices)
                np.save(join(output_directory, classification_str + '__intersect'), intersect)

            if save_original_weights or save_features_image:
                weights_orig = features_weights(current_subjects, dual_coefficients[0], sv_indices)

            if save_original_weights:
                np.save(join(output_directory, classification_str + '__weights'), weights_orig)

            if save_features_image:
                output_image=weights_to_nifti(input_image_atlas,weights_orig)
                output_image.to_filename(join(output_directory,classification_str+'__weights.nii'))

            if save_subject_classification:
                save_subjects_prediction(current_subjects, current_diagnosis, y, y_hat, join(output_directory, classification_str + '__subjects.csv'))

            results[(dx1, dx2)] = evaluate_prediction(y, y_hat)

    results_to_csv(results, dx_filter, join(output_directory, 'resume' + ('_balanced' if balanced else '_not_balanced') + '.csv'))