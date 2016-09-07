# -*- coding: utf-8 -*-
"""
Created on Mon Jun  6 10:00:08 2016

@author: jacquemont
"""

import clinica.pipeline.postprocessing.NBS as NBS

def create_network_based_statistic_pipeline(list_of_connectome_1, list_of_connectome_2, test, nb_permutation, size, output_prefix, 
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
                                                              function=NBS.Network_Based_Statistics), name='network_based_statistics')
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