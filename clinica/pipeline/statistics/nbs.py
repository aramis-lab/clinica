# -*- coding: utf-8 -*-
"""
Created on Thu Jun  2 11:15:38 2016

@author: jacquemont
"""

def network_based_statistic_pipeline(list_of_connectome_1, list_of_connectome_2, test, nb_permutation, size, output_prefix,
                                            working_directory, datasink_directory, threshold=0.05, significancy=0.05, tail=1, save_all=False):

    """
    This function create a pipeline which performs Network Based Statistic in the way of Zalesky (Zalesky et al. 2012).

    INPUTNODE:
        list_of_connectome_1 (list): list of connectome path from the first group
        list_of_connectome_2 (list): list of connectome path from the second group

    PARAMETER:
        test (str) : statistical test to use. Sould be 'permutation', 't test' or 'mann-whitney'.
        nb_permutation (int): number of permutation to perform.
        size (str): Size type to compute for the modules. Should be 'extent', 'intensity' or 'log'.
        output_prefix (str): A prefix that is prepended to all output files
        threshold (float): Primary threshold to use to define the module.
        significancy (float): Secondary threshold to consider the module size significant.
        tail (int): Tail of the statistic test. 2 is recommanded (default value).
        save_all (bool): If set True, all the modules and their respective p_value are keep even if their are not significant.
            Default value is False.

    OUTPUTNODE:
        list_of_graph_matrix_path (list): List of module(s) individual connection matrix path
        module_size_p_values (list): List of module(s)' size P value (in the same order as list_of_graph_matrix)
        p_value_matrix (str): Matrix of all connections who reach the threshold. Values in the matrix correspond
            to 1-p_value
    """

    import nipype.interfaces.io as nio
    import nipype.interfaces.utility as niu
    import nipype.pipeline.engine as pe
    import os.path as op
    import clinica.pipeline.statistics.nbs as nbs


    # Inputs checking

    if len(list_of_connectome_1)<2 or len(list_of_connectome_2)<2:
        raise IOError('List of connectomes should be greater than 1')

    for connectome in list_of_connectome_1:
        if not op.exists(connectome):
            raise IOError('file {} does not exist'.format(connectome))

    for connectome in list_of_connectome_2:
        if not op.exists(connectome):
            raise IOError('file {} does not exist'.format(connectome))

    if test not in ['permutation','t test', 'mann-whitney']:
        raise IOError('Set test parametter to "permutation","t test", or "mann-whitney".')

    if type(nb_permutation)!=int:
        raise IOError('nb_permutation shoul be a int.')

    if size not in ['extent', 'intensity', 'log']:
        raise IOError('Set size parametter to "extent", "intensity" or "log".')

    if type(threshold)!=float:
        raise IOError('threshold shoul be a float.')
    else:
        if threshold>1 or threshold<0.:
            raise IOError('threshold shoul be between 0 and 1.')

    if type(significancy)!=float:
        raise IOError('significancy shoul be a float.')
    else:
        if significancy>1 or significancy<0.:
            raise IOError('significancy shoul be between 0 and 1.')

    if tail!=1 and tail!=2:
        raise IOError('Set tail parametter to 1 or 2.')

    if type(save_all)!=bool:
        raise IOError('Set save_all  parametter to False or True.')

    # Nodes definition

    inputnode = pe.Node(niu.IdentityInterface(fields=['list_of_connectome_1', 'list_of_connectome_2']),
                        name='inputnode')
    inputnode.inputs.list_of_connectome_1 = list_of_connectome_1
    inputnode.inputs.list_of_connectome_2 = list_of_connectome_2

    network_based_statistics = pe.Node(interface=niu.Function(input_names=['list_of_connectome_1', 'list_of_connectome_2', 'test',
                                                                           'nb_permutation', 'size', 'output_prefix', 'threshold',
                                                                           'significancy', 'tail', 'save_all'],
                                                              output_names=['list_of_graph_matrix_path', 'module_size_p_values', 'p_value_matrix'],
                                                              function=nbs.network_based_statistics), name='network_based_statistics')
    network_based_statistics.inputs.test = test
    network_based_statistics.inputs.nb_permutation = nb_permutation
    network_based_statistics.inputs.size = size
    network_based_statistics.inputs.output_prefix = output_prefix
    network_based_statistics.inputs.threshold = threshold
    network_based_statistics.inputs.significancy = significancy
    network_based_statistics.inputs.tail = tail
    network_based_statistics.inputs.save_all = save_all

    outputnode = pe.Node(niu.IdentityInterface(fields=['list_of_graph_matrix_path', 'module_size_p_values', 'p_value_matrix']), name='outputnode')

    datasink = pe.Node(nio.DataSink(), name='datasink')
    datasink.inputs.base_directory = op.join(datasink_directory,'network_based_statistics/')

    # Build workflow

    wf = pe.Workflow(name='NBS')
    wf.base_dir = working_directory

    wf.connect([(inputnode, network_based_statistics, [('list_of_connectome_1', 'list_of_connectome_1'),
                                                       ('list_of_connectome_2','list_of_connectome_2')])])
    wf.connect([(network_based_statistics, outputnode, [('list_of_graph_matrix_path','list_of_graph_matrix_path'),
                                                        ('module_size_p_values','module_size_p_values'),
                                                        ('p_value_matrix','p_value_matrix')])])
    wf.connect([(network_based_statistics, datasink, [('list_of_graph_matrix_path','list_of_graph_matrix_path'),
                                                        ('module_size_p_values','module_size_p_values'),
                                                        ('p_value_matrix','p_value_matrix')])])

    return wf


def network_based_statistics(list_of_connectome_1, list_of_connectome_2, test, nb_permutation, size, output_prefix,
                             threshold=0.05, significancy=0.05, tail=1, save_all=False):
    """
    This function performs Network Based Statistic in the way of Zalesky (Zalesky et al. 2012).
    This is a Python implementation of Zalesky NBS. The statistic test used here can be one or two tail(s) and can be a permutation test, a T test or a
    Mann-Withney U test.
    take into account the existance of zeros in structural connexion matrix.
    If tail=1 is choose this function performs the contrast list_of_connectome_1 > list_of_connectome_2

    Args:
        list_of_connectome_1 (list): list of connectome path from the first group
        list_of_connectome_2 (list): list of connectome path from the second group
        test (str) : statistical test to use. Sould be 'permutation', 't test' or 'mann-whitney'.
        nb_permutation (int): number of permutation to perform.
        size (str): Size type to compute for the modules. Should be 'extent', 'intensity' or 'log'.
        output_prefix (str): A prefix that is prepended to all output files
        threshold (float): Primary threshold to use to define the module.
        significancy (float): Secondary threshold to consider the module size significant.
        tail (int): Tail of the statistic test. 2 is recommanded (default value).
        save_all (bool): If set True, all the modules and their respective p_value are keep even if their are not significant.
            Default value is False.

    Returns:
        list_of_graph_matrix_path (list): List of module(s) individual connection matrix path
        module_size_p_values (list): List of module(s)' size P value (in the same order as list_of_graph_matrix)
        p_value_matrix (str): Matrix of all connections who reach the threshold. Values in the matrix correspond
            to 1-p_value
    """
    from clinica.utils.statistics import permutation_test, t_test, Mann_Whitney
    from clinica.utils.postprocessing.NBS_size import log_size, compute_graph_size
    import copy
    import numpy as np
    import networkx as nx
    import os

    # Inputs checking

    if len(list_of_connectome_1)<2 or len(list_of_connectome_2)<2:
        raise IOError('List of connectomes should be greater than 1')

    if test not in ['permutation','t test', 'mann-whitney']:
        raise IOError('Set test parametter to "permutation","t test", or "mann-whitney".')

    if type(nb_permutation)!=int:
        raise IOError('nb_permutation shoul be a int.')

    if type(output_prefix)!=str:
        raise IOError('nb_permutation shoul be a str.')

    if size not in ['extent', 'intensity', 'log']:
        raise IOError('Set size parametter to "extent", "intensity" or "log".')

    if type(threshold)!=float:
        raise IOError('threshold shoul be a float.')
    else:
        if threshold>1 or threshold<0.:
            raise IOError('threshold shoul be between 0 and 1.')

    if type(significancy)!=float:
        raise IOError('significancy shoul be a float.')
    else:
        if significancy>1 or significancy<0.:
            raise IOError('significancy shoul be between 0 and 1.')

    if tail!=1 and tail!=2:
        raise IOError('Set tail parametter to 1 or 2.')

    if type(save_all)!=bool:
        raise IOError('Set save_all  parametter to False or True.')

    # initializing variables

    nb_subjects_1 = len(list_of_connectome_1)
    nb_subjects_2 = len(list_of_connectome_2)
    Concatenate_connectome = np.recfromtxt(list_of_connectome_1[0], delimiter=';', dtype=float)
    nb_connexion = (Concatenate_connectome.shape[0]*(Concatenate_connectome.shape[0]-1))/2
    connexion_index_vector = []
    vectorized_data = np.zeros((nb_subjects_1+nb_subjects_2, nb_connexion))
    permutation_max_size_vectors = np.zeros((nb_permutation))

    # Concatenation of all the connectome from the two list in one 3D array in the list order

    for connectome in list_of_connectome_1[1:]:
        Connect = np.recfromtxt(connectome, delimiter=';', dtype=float)
        Concatenate_connectome = np.dstack((Concatenate_connectome, Connect))
    for connectome in list_of_connectome_2:
        Connect = np.recfromtxt(connectome, delimiter=';', dtype=float)
        Concatenate_connectome = np.dstack((Concatenate_connectome, Connect))

    # Matrix for each subject from Concatenate_connectome is vectorized. The list connexion_index_vector allows to
    # come back to the matrix format. Her n is the number ID of the connection (i.e. the nth connection in vectorized
    # have for coordinate in the matrix format connexion_index_vector[n] )

    for k in np.arange(0, Concatenate_connectome.shape[2]):
        n=0
        if k==0:
            for i in np.arange(0, Concatenate_connectome.shape[0]):
                for j in np.arange(0,i):
                    connexion_index_vector += [(j,i)]
        for i in np.arange(0, Concatenate_connectome.shape[0]):
            for j in np.arange(0,i):
                vectorized_data[k,n] = Concatenate_connectome[j,i,k]
                n+=1

    # doing permutations of the data and permutation test/T Test on the permuted data on the fly.

    p_values = np.ones((1))
    graph_size = np.ones((1))
    tmp_count = 0.

    for i in range(nb_permutation+1):
        if i==0:
            data_perm = vectorized_data
        else:
            vector_perm= np.random.permutation(nb_subjects_1+nb_subjects_2).tolist()
            data_perm = vectorized_data[vector_perm,:]
        p_values_matrix = np.zeros(Concatenate_connectome[:,:,0].shape)
        for k in range(nb_connexion):
            list_1 = data_perm[:nb_subjects_1, k]
            list_2 = data_perm[nb_subjects_1:, k]
            if test=='permutation':
                p_value = permutation_test(list_1,list_2, nb_permutation, tails=tail)
            elif test=='t test':
                p_value = t_test(list_1,list_2, tails=tail)
            else:
                p_value = Mann_Whitney(list_1,list_2, tails=tail)
            if p_value<threshold:
                p_values_matrix[connexion_index_vector[k][0],connexion_index_vector[k][1]] = 1-p_value
        if i==0:
                G = nx.Graph(p_values_matrix)
                p_value_matrix_data = copy.deepcopy(p_values_matrix)
                list_of_sub_graph = list(nx.connected_component_subgraphs(G))
                p_values = np.ones(len(list_of_sub_graph))
                graph_size = np.zeros(len(list_of_sub_graph))
                tmp_count = 0.
                for nb, graph in enumerate(list_of_sub_graph):
                    graph_size[nb] = compute_graph_size(graph, size)
                print('|Permutation|      Max Size     |      Max Size     | Lowest  |')
                print('|           |       Random      |       Actual      | p-value |')
        else:
                maximum_size = 0
                G_rand = nx.Graph(p_values_matrix)
                for graph in list(nx.connected_component_subgraphs(G_rand)):
                    sub_size = compute_graph_size(graph, size)
                    if sub_size>maximum_size:
                        maximum_size=sub_size
                permutation_max_size_vectors[i-1] = maximum_size
                tmp_count += graph_size.max() < maximum_size
                tmp_p_value = tmp_count/float(i)
                print('|  {}/{}  |   {}   |   {}   |   {}   |'.format(i, nb_permutation, maximum_size, graph_size.max(), tmp_p_value))

    # Check for each component if it reach the significancy threshold and if yes, the subgraph adjacency matrix is stored

    for i, observed_size in enumerate(graph_size):
        k = 0.
        for random_size in permutation_max_size_vectors:
            k += observed_size < random_size
        p_values[i] = k/float(nb_permutation)

    list_of_graph_matrix = []
    list_of_p_value = []

    for i, p_val in enumerate(p_values):
        if save_all:
            sub_graph = nx.Graph(list_of_sub_graph[i])
            Adj = copy.deepcopy(p_value_matrix_data)
            for node in np.arange(0,p_value_matrix_data.shape[0]):
                 if node not in sub_graph.nodes():
                     Adj[node,:] = np.zeros(Adj[node,:].shape)
                     Adj[:,node] = np.zeros(Adj[:,node].shape)
            Adj = Adj.T + Adj
            list_of_graph_matrix += [Adj]
            list_of_p_value += [p_val]
        else:
            if p_val<significancy:
                sub_graph = nx.Graph(list_of_sub_graph[i])
                Adj = copy.deepcopy(p_value_matrix_data)
                for node in np.arange(0,p_value_matrix_data.shape[0]):
                    if node not in sub_graph.nodes():
                        Adj[node,:] = np.zeros(Adj[node,:].shape)
                        Adj[:,node] = np.zeros(Adj[:,node].shape)
                        Adj = Adj.T + Adj
                list_of_graph_matrix += [Adj]
                list_of_p_value += [p_val]

    working_directory = os.getcwd()

    os.mkdir(working_directory + '/' + output_prefix + '_Modules_matrix')

    list_of_graph_matrix_path = []

    for nb, module in enumerate(list_of_graph_matrix):
        new_path = working_directory + '/' + output_prefix + '_Modules_matrix/Module_' + str(nb) + '.csv'
        np.savetxt(new_path, module, delimiter=' ')
        list_of_graph_matrix_path += [new_path]

    module_size_p_values = working_directory + '/' + output_prefix + '_modules_size_p_values.csv'
    p_value_matrix = working_directory + '/' + output_prefix + '_matrix_p_values.csv'

    np.savetxt(module_size_p_values, np.array(list_of_p_value), delimiter=' ')
    np.savetxt(p_value_matrix, np.array(p_value_matrix_data), delimiter=' ')

    return list_of_graph_matrix_path, module_size_p_values, p_value_matrix
