# -*- coding: utf-8 -*-
"""
Created on Fri Apr 15 11:56:12 2016

@author: jacquemont
"""

def create_DTI_atlas_scalar_analysis(in_scalar_image, atlas_labels, atlas_scalar_image, working_directory, datasink_directory):
    
    """
    Perform tracts analysis according to a white matter atlas using a tensor-derived scalar image.

    This function perform the analysis of tracts using a white matter atlas and compute mean value of 
    the scalar on each tracts of this atlas. The function first coregister the subject scalar image 
    on the equivalent scalar image of the atlas and then use the labels to computs the statistics of the 
    scalar on each tracks of the whie matter atlas.

    Args:
        in_scalar_image (str): 3D image of the scalar obtained from the tensor
        atlas_labels (str): 3D Image of the white matter labels from the atlas 
        atlas_scalar_image (str): 3D image of the same scalar as in "in_scalar_image" but from the atlas

    Returns:
        out_stats_file (str): File containing for each tracts, the mean value of the scalar, the standart deviation
            and the nb of voxels.
    """    
    
    
    import nipype.interfaces.io as nio
    import nipype.interfaces.utility as niu
    import nipype.pipeline.engine as pe
    import os.path as op
    from clinica.pipeline.registration.mri_registration import antsRegistrationSyNQuick
            
    def DTI_atlas_scalar_analysis(input_image, atlas_labels_image):
    
        import nibabel as nib
        import numpy as np
        import pandas as pd
        import os.path as op
        
        outfile = op.abspath('scalar_stats.csv')
        
        DTI_atlas = nib.load(atlas_labels_image)
        atlas_image_data = DTI_atlas.get_data()
        
        labels = list(set(atlas_image_data.ravel()))
        stats_scalar = np.zeros((len(labels),4))
        
        In_image = nib.load(input_image)
        scalar_image_data = In_image.get_data()
        
        for index, label in enumerate(labels):
            stats_scalar[index, 0] = label
            atlas_label_index = np.array(np.where(atlas_image_data==label))
            nb_voxel = atlas_label_index.shape[1]
            stats_scalar[index, 3] = nb_voxel        
            labeled_voxel = labeled_voxel = scalar_image_data[atlas_label_index[0,:], atlas_label_index[1,:], atlas_label_index[2,:]]
            mean_scalar = labeled_voxel.mean()
            stats_scalar[index, 1] = mean_scalar
            std_scalar = labeled_voxel.std()
            stats_scalar[index, 2] = std_scalar
        
        columns = np.array(['Label', 'Mean scalar', 'Std scalar', 'Nb of voxel'])
        Data = pd.DataFrame(stats_scalar, columns=columns)
        
        Data.to_csv(outfile, sep=',', index=False)
    
        return outfile
        
    # Inputs existence checking

    inputs=[in_scalar_image, atlas_labels, atlas_scalar_image, working_directory, datasink_directory]     
        
    for input_file in inputs:
        if not op.exists(input_file):
            raise IOError('file {} does not exist'.format(input_file))
    
    # Nodes definition

    datasource = pe.Node(interface=nio.DataGrabber(infields=[], outfields=['in_scalar_image', 'atlas_labels', 'atlas_scalar_image']), name='datasource')
    datasource.inputs.template = '*'
    datasource.inputs.field_template = dict(in_scalar_image= in_scalar_image,
                                            atlas_labels=atlas_labels,
                                            atlas_scalar_image=atlas_scalar_image)
    datasource.inputs.template_args = dict(in_scalar_image=[[]],
                                           atlas_labels=[[]],
                                           atlas_scalar_image=[[]])
    datasource.inputs.sort_filelist = True
    
        
    inputnode = pe.Node(niu.IdentityInterface(fields=['in_scalar_image', 'atlas_labels', 'atlas_scalar_image']),
                        name='inputnode')
    
    antsRegistrationSyNQuick = pe.Node(interface=niu.Function(input_names=['fixe_image', 'moving_image'], output_names=['image_warped', 'affine_matrix', 'warp', 'inverse_warped', 'inverse_warp'],
                                                              function=antsRegistrationSyNQuick), name='antsRegistrationSyNQuick')
    
    scalar_analysis = pe.Node(interface=niu.Function(input_names=['input_image', 'atlas_labels_image'], output_names=['outfile'],
                                                              function=DTI_atlas_scalar_analysis), name='scalar_analysis')
    
    outputnode = pe.Node(niu.IdentityInterface(fields=['out_stats_file', 'affine_matrix', 'warp', 'inverse_warp', 'inverse_warped']), name='outputnode')

    datasink = pe.Node(nio.DataSink(), name='datasink')
    datasink.inputs.base_directory = op.join(datasink_directory,'dti_scalar_analysis/')
    
    # Building workflow
    
    wf = pe.Workflow(name='dti_scalar_analysis')
    wf.base_dir = working_directory
    
    wf.connect([(datasource, inputnode, [('in_scalar_image', 'in_scalar_image'),('atlas_labels', 'atlas_labels'),('atlas_scalar_image','atlas_scalar_image')])])
    wf.connect([(inputnode, antsRegistrationSyNQuick, [('in_scalar_image', 'moving_image'),('atlas_scalar_image', 'fixe_image')])])
    wf.connect([(inputnode, scalar_analysis, [('atlas_labels', 'atlas_labels_image')])])
    wf.connect([(antsRegistrationSyNQuick, scalar_analysis, [('image_warped', 'input_image')])])
    wf.connect([(antsRegistrationSyNQuick, outputnode, [('affine_matrix', 'affine_matrix'),('warp','warp'),('inverse_warp','inverse_warp'),('inverse_warped','inverse_warped')])])
    wf.connect([(antsRegistrationSyNQuick, datasink, [('affine_matrix', 'affine_matrix'),('warp','warp'),('inverse_warp','inverse_warp'),('inverse_warped','inverse_warped')])])
    wf.connect([(scalar_analysis, outputnode, [('outfile', 'out_stats_file')])])
    wf.connect([(scalar_analysis, datasink, [('outfile', 'out_stats_file')])])
    
    return wf
    