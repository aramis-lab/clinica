# -*- coding: utf-8 -*-
"""
Created on Mon Jun  6 10:00:08 2016

@author: jacquemont
"""

import clinica.pipeline.postprocessing.connection_wise_analysis as CWA

def create_connection_wise_analysis_pipeline(list_of_connectome_1, list_of_connectome_2, test, working_directory,
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
                                                              function=CWA.connection_wise_analysis), name='connection_wise_analysis')
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
