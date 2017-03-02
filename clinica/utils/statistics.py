# -*- coding: utf-8 -*-
"""
Created on Thu Jun  2 11:15:38 2016

@author: jacquemont
"""

def permutation_test(vector_1, vector_2, number_of_permutation, tails=2):
    """
    This function take two vectors as entry and return the p_value of a permutation test. 
    If permutation not done, p_value is set to 1.
    This function performs a two tails permutation on *non-nul values*
    If tails is set to 1 the contrast is supposed to be vector_1 > vector_2
    
    WARNING: here 0 are not removed from the dataset
    If tails is set to 1 the contrast is supposed to be vector_1 > vector_2
    """
    
    import numpy as np
    
    def exact_mc_perm_test(xs, ys, nmc, tails=2):
        n, k = len(xs), 0.
        if tails==2:
            diff = np.abs(np.mean(xs) - np.mean(ys))
        if tails==1:
            diff = np.mean(xs) - np.mean(ys)
        zs = np.concatenate([xs, ys])
        for j in range(nmc):
            np.random.shuffle(zs)
            if tails==2:
                diff_random = np.abs(np.mean(zs[:n]) - np.mean(zs[n:]))
            if tails==1:
                diff_random = np.mean(zs[:n]) - np.mean(zs[n:])
            k += diff < diff_random
        return k/float(nmc)
    
    p_value = 1.
    
    X = vector_1
    Y = vector_2
    
    if X!=np.array([]) and Y!=np.array([]):   
        if X.sum()!=0 or Y.sum()!=0:
            p_value = exact_mc_perm_test(X,Y,number_of_permutation, tails)
    
    return p_value

def t_test(vector_1, vector_2, tails=2):
    """
    This function take two vectors as entry and return the p_value. of a t_test. 
    If the t_test not done (only zeros in both vectors), p_value is set to 1.
    This function performs a t_test on *non-nul values*.
    If tails is set to 1 the contrast is supposed to be vector_1 > vector_2.
    
    WARNING: here 0 are not removed from the dataset
    If tails is set to 1 the contrast is supposed to be vector_1 > vector_2
    """
    
    import numpy as np
    from scipy.stats import ttest_ind
    
    p_value = 1.
    
    X = vector_1
    Y = vector_2
    
    if X!=np.array([]) and Y!=np.array([]):
        p_value= ttest_ind(X,Y,equal_var=False).pvalue
        if tails==1:
            array_X = np.array(X)
            array_Y = np.array(Y)
            if array_X.mean()>array_Y.mean():
                p_value = float(p_value)/2.0
            else:
                p_value = 1 - float(p_value)/2.0
    
    return p_value
    
def Mann_Whitney(vector_1, vector_2, tails=2):
    """
    This function take two vectors as entry and return 
    the p_value of a  Mann_Whitney U test. 
    If the U test not done (only 0 in both vectors), p_value is set to 1.
    This function performs a U test on *non-nul values*
    
    WARNING: here 0 are not removed from the dataset
    If tails is set to 1 the contrast is supposed to be vector_1 > vector_2
    """
    
    import numpy as np
    from scipy.stats import mannwhitneyu
    
    p_value = 1.
    
    X = vector_1
    X_non_nul = [X[k] for k in np.nonzero(X)[0].tolist()]
 
    Y = vector_2
    Y_non_nul = [Y[k] for k in np.nonzero(Y)[0].tolist()]
    
    if X!=np.array([]) and Y!=np.array([]):
        if X_non_nul!=np.array([]) or Y_non_nul!=np.array([]):
            p_value= mannwhitneyu(X,Y).pvalue
            if tails==2:
                p_value = 2*p_value
    
    return p_value
    
def fdr_correction_matrix(p_value_matrix, template=None):
    """
    This function take a p value matrix as entry and return the p_value corrected for False Rate Discovery. 
    If not all statistical tests have been performed (typically in DTI at a absent connection) a template
    matrix (which is a binary metrix with 1 if the test is performed and 0 else) with the same shape as 
    p_value_matrix input of the actually performed test can be provide at input.
    """
    
    import numpy as np
    from mne.stats import fdr_correction
    
    if type(template)==type(p_value_matrix):
        if p_value_matrix.shape!=template.shape:
            raise IOError('p_value_matrix and template should have the same shape.')
         
    if type(template)==type(p_value_matrix):  
        p_value_corrected = np.ones(p_value_matrix.shape)
        reject_test = np.zeros(p_value_matrix.shape, dtype=bool)     
        eff_p_value = []
        index_of_eff_p_value = []    
        for i in np.arange(0, p_value_matrix.shape[0]):
            for j in np.arange(0,i):
                if template[j,i]==1:
                    eff_p_value += [p_value_matrix[j,i]]
                    index_of_eff_p_value += [(j,i)]     
        reject, p_corrected = fdr_correction(eff_p_value)     
        for i, corrected in enumerate(p_corrected):
            p_value_corrected[index_of_eff_p_value[i][0],index_of_eff_p_value[i][1]] = corrected
            reject_test[index_of_eff_p_value[i][0],index_of_eff_p_value[i][1]] = reject[i]
    elif not template:
         reject_test, p_value_corrected = fdr_correction(p_value_matrix)
    else:
         raise IOError('template input should be an numpy array or None.') 
    return reject_test, p_value_corrected

def create_new_feature_tsv(subjects_visits_tsv, bids_dir, dest_tsv, added_features):
    """This func is to add new features(columns) from the subjects_visits_list tsv file, and use the generated file in the statistical analysis"""
    import pandas as pd
    from os.path import join, isfile
    import logging
    from pandas import concat
    logging.basicConfig(level=logging.DEBUG)

    if not isfile(join(bids_dir, 'participants.tsv')):
        raise Exception('participants.tsv not found')
    sub_set = pd.io.parsers.read_csv(subjects_visits_tsv, sep='\t')
    all_set = pd.io.parsers.read_csv(join(bids_dir, 'participants.tsv'), sep='\t')
    selected_subj = all_set[all_set.participant_id.isin(list(sub_set.participant_id))]
    if len(sub_set.participant_id) != len(selected_subj.participant_id):
        missing_subj = set(list(sub_set.participant_id)) - set(list(selected_subj.participant_id))
        msg = "The missing subjects are %s" % list(missing_subj)
        logging.info(msg)
        raise Exception('There are subjects which are not included in full dataset, please check it out')

    new_features = selected_subj[added_features]
    all_features = concat([sub_set, new_features], axis=1)
    all_features.to_csv(dest_tsv, sep='\t', index=False)