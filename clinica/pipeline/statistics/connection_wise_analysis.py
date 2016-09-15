# -*- coding: utf-8 -*-
"""
Created on Thu Jun  2 11:15:38 2016

@author: jacquemont
"""

def connection_wise_analysis_pipeline(list_of_connectome_1, list_of_connectome_2, test, working_directory,
                                            datasink_directory, FDR_correction=True, tail=1, nb_permutation=0):

    """
    This function create a pipeline which performs connection wise analysis.
    The statistic test used here can be one or two tail(s) and can be a permutation test, a T test or a
    Mann-Withney U test.
    If tail=1 is choose this function performs the contrast list_of_connectome_1 > list_of_connectome_2

    INPUTNODE:
        list_of_connectome_1 (list): list of connectome path from the first group
        list_of_connectome_2 (list): list of connectome path from the second group

    PARAMETER:
        test (str) : statistical test to use. Sould be 'permutation', 't test' or 'mann-whitney'.
        FDR_correction (bool): Perform an FDR correction if True. Default value is True
        tail (int): Tail of the statistic test. 2 is recommanded (default value).
        nb_permutation (int): number of permutation to perform.

    OUTPUTNODE:
        p_values_path (str): Uncorrected P value matrix
        p_values_corrected_path(str - optional): FDR corrected P value matrix (if FDR_correction set to True).

    """

    import nipype.interfaces.io as nio
    import nipype.interfaces.utility as niu
    import nipype.pipeline.engine as pe
    import os.path as op

    import clinica.pipeline.statistics.connection_wise_analysis as cwa


    # Inputs checking

    if len(list_of_connectome_1)<2 or len(list_of_connectome_2)<2:
        raise IOError('List of connectomes should be greater than 1')

    for connectome in list_of_connectome_1:
        if not op.exists(connectome):
            raise IOError('file {} does not exist'.format(connectome))

    for connectome in list_of_connectome_2:
        if not op.exists(connectome):
            raise IOError('file {} does not exist'.format(connectome))

    if type(nb_permutation)!=int:
        raise IOError('nb_permutation shoul be a int.')

    if test not in ['permutation','t test', 'mann-whitney']:
        raise IOError('Set test parametter to "permutation","t test", or "mann-whitney".')
    elif test=='permutation':
        if nb_permutation<1:
            raise IOError('nb_permutation shoul be > 0.')

    if tail!=1 and tail!=2:
        raise IOError('Set tail parametter to 1 or 2.')

    # Nodes definition

    inputnode = pe.Node(niu.IdentityInterface(fields=['list_of_connectome_1', 'list_of_connectome_2']),
                        name='inputnode')
    inputnode.inputs.list_of_connectome_1 = list_of_connectome_1
    inputnode.inputs.list_of_connectome_2 = list_of_connectome_2

    connection_wise_analysis = pe.Node(interface=niu.Function(input_names=['list_of_connectome_1', 'list_of_connectome_2', 'test',
                                                                           'FDR_correction', 'tail', 'nb_permutation'],
                                                              output_names=['p_values_path', 'p_values_corrected_path'],
                                                              function=cwa.connection_wise_analysis), name='connection_wise_analysis')
    connection_wise_analysis.inputs.test = test
    connection_wise_analysis.inputs.FDR_correction = FDR_correction
    connection_wise_analysis.inputs.nb_permutation = nb_permutation
    connection_wise_analysis.inputs.tail = tail

    outputnode = pe.Node(niu.IdentityInterface(fields=['p_values_path', 'p_values_corrected_path']), name='outputnode')

    datasink = pe.Node(nio.DataSink(), name='datasink')
    datasink.inputs.base_directory = op.join(datasink_directory,'connection_wise_analysis/')

    # Build workflow

    wf = pe.Workflow(name='Connection_wise_analysis')
    wf.base_dir = working_directory

    wf.connect([(inputnode, connection_wise_analysis, [('list_of_connectome_1', 'list_of_connectome_1'),
                                                       ('list_of_connectome_2','list_of_connectome_2')])])
    wf.connect([(connection_wise_analysis, outputnode, [('p_values_path','p_values_path'),
                                                        ('p_values_corrected_path','p_values_corrected_path')])])
    wf.connect([(connection_wise_analysis, datasink, [('p_values_path','p_values_path')])])

    if FDR_correction:
        wf.connect([(connection_wise_analysis, datasink, [('p_values_corrected_path','p_values_corrected_path')])])

    return wf


def connection_wise_analysis(list_of_connectome_1, list_of_connectome_2, test, FDR_correction, tail=1, nb_permutation=0):
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
    from clinica.pipeline.postprocessing.connectome_stats import permutation_test, t_test, Mann_Whitney, fdr_correction_matrix
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

    Concatenate_connectome_1 = np.recfromtxt(list_of_connectome_1[0], delimiter=' ', dtype=float)

    if len(Concatenate_connectome_1)<2:
        raise IOError('Connectome matrix csv from first list file should be separated by a space ' '.')
    if len(Concatenate_connectome_1)<2:
        raise IOError('Connectome matrix csv from second list file should be separated by a space ' '.')

    Concatenate_connectome_2 = np.recfromtxt(list_of_connectome_2[0], delimiter=' ', dtype=float)

    # Concatenation of all the connectome from the two list in two 3D array in the list order

    for connectome in list_of_connectome_1[1:]:
        Connect = np.recfromtxt(connectome, delimiter=' ', dtype=float)
        Concatenate_connectome_1 = np.dstack((Concatenate_connectome_1, Connect))

    for connectome in list_of_connectome_2[1:]:
        Connect = np.recfromtxt(connectome, delimiter=' ', dtype=float)
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
                    p_value = permutation_test(list_1,list_2, nb_permutation, tails=tail)
                elif test=='t test':
                    p_value = t_test(list_1,list_2, tails=tail)
                else:
                    p_value = Mann_Whitney(list_1,list_2, tails=tail)
                p_values[j,i] = p_value

    if FDR_correction:
        reject, p_values_corrected = fdr_correction_matrix(p_values, template=effected_test)

    working_directory = os.getcwd()

    p_values_path = working_directory + '/P_values_matrix.csv'
    np.savetxt(p_values_path, p_values, delimiter=' ')

    if FDR_correction:
        p_values_corrected_path = working_directory + '/P_values_matrix_FDR_corrected.csv'
        np.savetxt(p_values_corrected_path, p_values_corrected, delimiter=' ')
    else:
        p_values_corrected_path = None

    return p_values_path, p_values_corrected_path
