# -*- coding: utf-8 -*-
"""
Created on Thu Jun  2 11:15:38 2016

@author: jacquemont
"""

import clinica.pipeline.postprocessing.connectome_stats as connectome_stats

def Connection_wise_analysis(list_of_connectome_1, list_of_connectome_2, test, FDR_correction, tail=1, nb_permutation=0):
    """
    This function performs connection wise analysis. 
    The statistic test used here can be one or two tail(s) and can be a permutation test, a T test or a
    Mann-Withney U test. 
    If tail=1 is choose this function performs the contrast list_of_connectome_1 > list_of_connectome_2
    
    Args:
        list_of_connectome_1 (list): list of connectome path from the first group
        list_of_connectome_2 (list): list of connectome path from the second group
        test (str) : statistical test to use. Sould be 'permutation', 't test' or 'mann-whitney'.
        FDR_correction (bool): Perform an FDR correction if True. Default value is True
        tail (int): Tail of the statistic test. 2 is recommanded (default value).
        nb_permutation (int): number of permutation to perform.
        
    Returns:
        p_values_path (str): Uncorrected P value matrix
        p_values_corrected_path(str - optional): FDR corrected P value matrix (if FDR_correction set to True).  
        
    """
    import numpy as np
    import os

    # Inputs checking

    if len(list_of_connectome_1)<2 or len(list_of_connectome_2)<2:
        raise IOError('List of connectomes should be greater than 1')
        
    if test not in ['permutation','t test', 'mann-whitney']:
        raise IOError('Set test parametter to "permutation","t test", or "mann-whitney".')
    
    if type(nb_permutation)!=int:
        raise IOError('nb_permutation shoul be a int.')
    elif test=="permutation" and nb_permutation<1:
        raise IOError('nb_permutation shoul superior to 0.')
    
    if type(FDR_correction)!=bool:
        raise IOError('Set FDR_correction parameter to True or False.')

    if tail!=1 and tail!=2:
        raise IOError('Set tail parametter to 1 or 2.')
    
        
    # initializing variables
        
    Concatenate_connectome_1 = np.recfromtxt(list_of_connectome_1[0], delimiter=';', dtype=float)
    Concatenate_connectome_2 = np.recfromtxt(list_of_connectome_2[0], delimiter=';', dtype=float)    
    
    # Concatenation of all the connectome from the two list in two 3D array in the list order  
    
    for connectome in list_of_connectome_1[1:]: 
        Connect = np.recfromtxt(connectome, delimiter=';', dtype=float)
        Concatenate_connectome_1 = np.dstack((Concatenate_connectome_1, Connect))
        
    for connectome in list_of_connectome_2[1:]: 
        Connect = np.recfromtxt(connectome, delimiter=';', dtype=float)
        Concatenate_connectome_2 = np.dstack((Concatenate_connectome_2, Connect))

    if Concatenate_connectome_1[:,:,0].shape!=Concatenate_connectome_1[:,:,0].shape:
        raise IOError('Connectomes from list 1 and list 2 should have the same shape.')
        
    # doing permutations of the data and permutation test/T Test on the permuted data on the fly. 
    
    p_values = np.ones(Concatenate_connectome_1[:,:,0].shape)
    effected_test = np.zeros(Concatenate_connectome_1[:,:,0].shape)
    
    for i in np.arange(0, Concatenate_connectome_1[0,:,0].shape[0]):
        for j in np.arange(0,i):
            list_1 = Concatenate_connectome_1[j,i,:]
            list_2 = Concatenate_connectome_2[j,i,:]
            if list_1.sum()!=0 or list_2.sum()!=0:
                effected_test[j,i] = 1
                if test=='permutation':
                    p_value = connectome_stats.permutation_test(list_1,list_2, nb_permutation, tails=tail)
                elif test=='t test':
                    p_value = connectome_stats.t_test(list_1,list_2, tails=tail)
                else:
                    p_value = connectome_stats.Mann_Whitney(list_1,list_2, tails=tail)
                p_values[j,i] = p_value
    
    if FDR_correction:
        reject, p_values_corrected = connectome_stats.fdr_correction_matrix(p_values, template=effected_test)
    
    working_directory = os.getcwd()    
    
    p_values_path = working_directory + '/P_values_matrix.csv'
    np.savetxt(p_values_path, p_values, delimiter=';')
    
    if FDR_correction:
        p_values_corrected_path = working_directory + '/P_values_matrix_FDR_corrected.csv'
        np.savetxt(p_values_corrected_path, p_values_corrected, delimiter=';')
    else:
        p_values_corrected_path = None
    
    return p_values_path, p_values_corrected_path
            
