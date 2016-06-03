# -*- coding: utf-8 -*-
"""
Created on Thu Jun  2 11:15:38 2016

@author: jacquemont
"""

def NBS_permutation_test(vector_1, vector_2, number_of_permutation, tails=2):
    """
    This function take two vectors as entry and return the p_value of a permutation test. 
    If permutation not done, p_value is set to 1.
    This function performs a two tails permutation on *non-nul values*
    If tails is set to 1 the contrast is supposed to be vector_1 > vector_2
    
    WARNING: here 0 are not removed from the dataset
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
            k += diff <= diff_random # EDIT: must be <= and not only < because if not, unexisting connexion (only 0) become significant
        return k/float(nmc)
    
    p_value = 1.
    
    X = vector_1
    Y = vector_2
    
    if X!=[] and Y!=[]:
        p_value = exact_mc_perm_test(X,Y,number_of_permutation, tails)
    
    return p_value

def NBS_t_test(vector_1, vector_2, tails=2):
    """
    This function take two vectors as entry and return the p_value. of a t_test. 
    If the t_test not done (only zeros in both vectors), p_value is set to 1.
    This function performs a t_test on *non-nul values*.
    If tails is set to 1 the contrast is supposed to be vector_1 > vector_2.
    
    WARNING: here 0 are not removed from the dataset
    """
    
    import numpy as np
    from scipy.stats import ttest_ind
    
    p_value = 1.
    
    X = vector_1
    Y = vector_2
    
    if X!=[] and Y!=[]:
        p_value= ttest_ind(X,Y,equal_var=False).pvalue
        if tails==1:
            array_X = np.array(X)
            array_Y = np.array(Y)
            if array_X.mean()>array_Y.mean():
                p_value = float(p_value)/2.0
            else:
                p_value = 1 - float(p_value)/2.0
    
    return p_value
    
def NBS_Mann_Whitney(vector_1, vector_2, tails=2):
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
    
    if X!=[] and Y!=[]:
        if X_non_nul!=[] or Y_non_nul!=[]:
            p_value= mannwhitneyu(X,Y).pvalue
            if tails==2:
                p_value = 2*p_value
    
    return p_value