# -*- coding: utf-8 -*-
"""
Created on Wed Jun  1 10:05:42 2016

@author: jacquemont
"""

def Connectome_construction_pipeline(in_parcellation, configuration_file, lut_type, lut_path, in_tracks, connectome_metric, working_directory, datasink_directory, in_scalar_image='', zeros_diagonal=True):
    
    """
    Perform the construction of the subject connectome.

    This function combine tracks obtained from tractography with a parcellated cortex in order to compute the 
    subject connectome as a connection matrix. This function first prepare a parcellated image for connectome 
    construction by modifying the image values; typically this involves making the parcellation intensities 
    increment from 1 to coincide with rows and columns of a matrix. The configuration file passed as the second 
    argument specifies the indices that should be assigned to different structures; examples of such configuration 
    files are provided in the MRtrix3 directory in src//dwi//tractography//connectomics//example_configs//. Then 
    a connectome matrix is generated from a streamlines file and the configurate node parcellation image.

    Args:
        in_parcellation (str): 3D image of the parcellated cortex
        configuration_file (str): the MRtrix connectome configuration file specifying desired nodes & indices 
        lut_type (str): Type of lookup table should be 'freesurfer', 'aal', 'itksnap' or 'basic' otherwise
        lut_path (str): Path of the lookup table. If type is 'FreeSurfer' then typically 'FreeSurferColorLUT.txt',
            if it is 'aal' then typically 'ROI_MNI_V4.txt', if it is 'itksnap' the ntypically this include the ITK atlas
            'LUT_GM.txt'. Otherwise, if it is 'basic' then this is a basic lookup table consisting of index / name pairs
        in_tracks (str): the input track file
        connectome_metric (str): specify the edge weight metric. Options are: count (default), meanlength, invlength, 
            invnodevolume, invlength_invnodevolume, mean_scalar. If mean_scalar is selected an associate scalar image
            is espected in the in_scalar_image
        in_scalar_image (str): provide the associated image for the mean_scalar metric
        zero_diagonal (bool): If set True set all diagonal entries in the matrix to zero (default is True)

    Returns:
        out_connectome (str): The output .csv file containing edge weights store in a matrix.
    """    
    
    import nipype.interfaces.io as nio
    import nipype.interfaces.utility as niu
    import nipype.pipeline.engine as pe
    import os.path as op
    
    def labelconfig(input_parcellation, config_in, lut_type, lut_path):

        import subprocess
        import os.path as op
        
        if lut_type.lower() in ['freesurfer', 'aal', 'itksnap']:
            if lut_type.lower()=='freesurfer':
                lut='-lut_freesurfer'
            if lut_type.lower()=='aal':
                lut='-lut_aal'
            if lut_type.lower()=='itksnap':
                lut='-lut_itksnap'
        else:
            lut='-lut_basic'

        out_image = op.abspath('parcellation_config.nii')

        cmd = 'labelconfig ' + lut + ' ' + lut_path + ' ' + input_parcellation + ' ' + config_in + ' ' + out_image
        subprocess.call([cmd], shell=True)

        return out_image
        
    def tck2connectome(in_parcellation, in_tracks, metric='count', scalar_image='', zero_diagonal=True):

        import subprocess
        import os.path as op

        out_connectome = op.abspath('connectome.csv')
        
        if zero_diagonal:
            diag = '-zero_diagonal '
        else:
            diag = ''
            
        if metric=='mean_scalar':
            if scalar_image:
                image = '-image ' + scalar_image + ' '
            else:
                raise IOError('A metric "mean_scalar" should be associated with a valid scalar_image input.')
        else:
            image=''

        cmd = ' tck2connectome ' + diag + ' -metric ' + metric + ' ' + image + in_tracks + ' ' + in_parcellation + ' ' + out_connectome
        subprocess.call([cmd], shell=True)

        if not op.exists(out_connectome):
            raise IOError('file {} does not exist, cmd is : {}'.format(out_connectome, cmd))
        return out_connectome

    # Inputs existence checking

    inputs=[in_parcellation, configuration_file, lut_path, in_tracks, working_directory, datasink_directory]     
        
    for input_file in inputs:
        if not op.exists(input_file):
            raise IOError('file {} does not exist'.format(input_file))
    
    if lut_type.lower() not in  ['freesurfer', 'aal', 'itksnap', 'basic']:
        raise IOError('lut_type should be in "freesurfer", "aal", "itksnap" or "basic"')
    
    if connectome_metric.lower() not in  ['count', 'meanlength', 'invlength', 'invnodevolume', 'invlength_invnodevolume', 'mean_scalar']:
        raise IOError('connectome_metric should be in "count", "meanlength", "invlength", "invnodevolume", "invlength_invnodevolume", "mean_scalar"')
    elif connectome_metric.lower()=='mean_scalar':
        if not op.exists(in_scalar_image):
            raise IOError('A "mean_scalar" connectome_metric value should be associated with a valid in_scalar_image input. file {} does not exist'.format(in_scalar_image))
            
    if type(zeros_diagonal)!=bool:
        raise IOError('zeros_diagonal should be set to "True" or "False".')
    
    # Nodes definition

    datasource = pe.Node(interface=nio.DataGrabber(infields=[], outfields=['in_parcellation', 'configuration_file','lut_path', 'in_tracks']), name='datasource')
    datasource.inputs.template = '*'
    datasource.inputs.field_template = dict(in_parcellation= in_parcellation,
                                            configuration_file=configuration_file,
                                            lut_path=lut_path,
                                            in_tracks=in_tracks)
    datasource.inputs.template_args = dict(in_parcellation=[[]],
                                           configuration_file=[[]],
                                           lut_path=[[]],
                                           in_tracks=[[]])
    datasource.inputs.sort_filelist = True
    
        
    inputnode = pe.Node(niu.IdentityInterface(fields=['in_parcellation', 'configuration_file', 'lut_path', 'in_tracks']),
                        name='inputnode')
    
    label_config = pe.Node(interface=niu.Function(input_names=['input_parcellation', 'config_in', 'lut_type', 'lut_path'], output_names=['out_image'],
                                                              function=labelconfig), name='label_config')
    label_config.inputs.lut_type = lut_type
    
    tck_to_connectome = pe.Node(interface=niu.Function(input_names=['in_parcellation', 'in_tracks', 'metric', 'scalar_image', 'zero_diagonal'], output_names=['out_connectome'],
                                                              function=tck2connectome), name='tck_to_connectome')
    tck_to_connectome.inputs.zeros_diagonal = zeros_diagonal
    tck_to_connectome.inputs.scalar_image = in_scalar_image
    tck_to_connectome.inputs.metric = connectome_metric
                                                              
    outputnode = pe.Node(niu.IdentityInterface(fields=['out_connectome']), name='outputnode')

    datasink = pe.Node(nio.DataSink(), name='datasink')
    datasink.inputs.base_directory = op.join(datasink_directory,'connectome/')
    
    wf = pe.Workflow(name='compute_connectome')
    wf.base_dir = working_directory
    
    wf.connect([(datasource, inputnode, [('in_parcellation','in_parcellation'), 
                                         ('configuration_file','configuration_file'),
                                         ('lut_path','lut_path'),
                                         ('in_tracks','in_tracks')])])
    wf.connect([(inputnode, label_config, [('in_parcellation','input_parcellation'),
                                           ('configuration_file','config_in'),
                                           ('lut_path','lut_path')])])
    wf.connect([(inputnode, tck_to_connectome, [('in_tracks','in_tracks')])])
    wf.connect([(label_config, tck_to_connectome, [('out_image','in_parcellation')])])
    wf.connect([(tck_to_connectome, outputnode, [('out_connectome','out_connectome')])])
    wf.connect([(tck_to_connectome, datasink, [('out_connectome','out_connectome')])])
    
    return wf