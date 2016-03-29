# -*- coding: utf-8 -*-
"""
Created on Thu Jan 28 16:57:18 2016

@author: jacquemont
"""

def tractography_pipeline(datasink_directory, tractography_algorithm='iFOD2', nb_track=150000, step_size_tracking=0.12, cutoff_value_tracking=0.1, n_harmonics=6, name='tractography' ):
    
    import nipype.interfaces.io as nio           # Data i/o
    import nipype.interfaces.utility as niu     # utility
    import nipype.pipeline.engine as pe          # pypeline engine
    import nipype.interfaces.mrtrix as mrtrix
    #import nipype.interfaces.mrtrix3 as mrtrix3
    import os
    import os.path as op
    
    def maskfilter(in_mask, filter_type, npass):
        
        import subprocess
        import os.path as op
        
        inputs=[in_mask]     
        
        for input_file in inputs:
            try:
                with open(input_file): pass
            except IOError:
                print "file {} does not exist".format(input_file)       
        
        out_file = op.abspath('out_mask.nii.gz')
        
        cmd = 'maskfilter ' + in_mask + ' ' + filter_type + ' out_mask.nii.gz -npass ' + str(npass)    
        subprocess.call([cmd], shell=True)      
        
        return out_file
    
    def dwi2tensor(in_dwi, in_bvecs, in_bvals, in_mask):
        
        import subprocess
        import os.path as op

        inputs=[in_dwi, in_bvecs, in_bvals, in_mask]     
        
        for input_file in inputs:
            try:
                with open(input_file): pass
            except IOError:
                print "file {} does not exist".format(input_file) 
                
        tensor_file = op.abspath('tensor.mif')
        
        cmd = 'dwi2tensor -mask ' + in_mask + ' -fslgrad ' + in_bvecs + ' ' +  in_bvals + ' ' + in_dwi + ' tensor.mif'    
        subprocess.call([cmd], shell=True)
    
        return tensor_file
        
    def tensor2fa(tensor_file, in_mask):
        
        import subprocess
        import os.path as op
        
        inputs=[tensor_file, in_mask]     
        
        for input_file in inputs:
            try:
                with open(input_file): pass
            except IOError:
                print "file {} does not exist".format(input_file)

        FA_map = op.abspath('fa.mif')
        
        cmd = ' tensor2metric -mask ' + in_mask + ' -fa fa.mif ' + tensor_file
        subprocess.call([cmd], shell=True)
    
        return FA_map
        
        
    def dwi2response(in_dwi, in_bvecs, in_bvals, in_mask, lmax=6):
        
        import subprocess
        import os.path as op

        inputs=[in_dwi, in_bvecs, in_bvals, in_mask]     
        
        for input_file in inputs:
            try:
                with open(input_file): pass
            except IOError:
                print "file {} does not exist".format(input_file) 
                
        response_file = op.abspath('response_file.txt')
        
        cmd = 'dwi2response -fslgrad ' + in_bvecs + ' ' + in_bvals + ' -lmax ' + str(lmax) + ' -mask ' + in_mask + ' ' +  in_dwi + ' response_file.txt'    
        subprocess.call([cmd], shell=True)
    
        return response_file
        
    def dwi2fod(in_dwi, in_bvecs, in_bvals, in_mask, response_file, lmax=6):
        
        import subprocess
        import os.path as op
        
        inputs=[in_dwi, in_bvecs, in_bvals, in_mask, response_file]     
        
        for input_file in inputs:
            try:
                with open(input_file): pass
            except IOError:
                print "file {} does not exist".format(input_file)

        SH_coeff_image = op.abspath('SH_coeff_image.nii.gz')
        
        cmd = 'dwi2fod -fslgrad ' + in_bvecs + ' ' + in_bvals + ' -lmax ' + str(lmax) + ' -mask ' + in_mask + ' ' +  in_dwi + ' ' + response_file + ' SH_coeff_image.nii.gz'
        subprocess.call([cmd], shell=True)
    
        return SH_coeff_image
        
    def tckgen(algorithm, step_size, number_tracks, in_bvecs, in_bvals, white_matter_mask, SH_coeff_image):
        
        import subprocess
        import os.path as op
        
        inputs=[in_bvecs, in_bvals, white_matter_mask, SH_coeff_image]     
        
        for input_file in inputs:
            try:
                with open(input_file): pass
            except IOError:
                print "file {} does not exist".format(input_file)

        out_tracks = op.abspath('out_tracks.tck')
        
        cmd = 'tckgen -algorithm ' + algorithm + ' -step ' + str(step_size) + ' -number ' + str(number_tracks) + ' -fslgrad ' + in_bvecs + ' ' + in_bvals + ' -seed_image ' + white_matter_mask + ' ' + SH_coeff_image + ' out_tracks.tck'

        subprocess.call([cmd], shell=True)
    
        return out_tracks
        
    # Tractography by MRtrix3

    inputnode = pe.Node(niu.IdentityInterface(fields=['in_dwi', 'in_bvecs', 'in_bvals', 'B0_mask', 'white_matter_mask']),
                        name='inputnode')
    
    # Computation of the FA_map
    
    dwi2tensor = pe.Node(interface=niu.Function(input_names=['in_dwi', 'in_bvecs', 'in_bvals', 'in_mask'], 
                                                  output_names=['tensor_file'], function=dwi2tensor), name='dwi2tensor')
                                                  
    tensor2fa_metric = pe.Node(interface=niu.Function(input_names=['tensor_file', 'in_mask'], 
                                                  output_names=['FA_map'], function=tensor2fa), name='tensor2fa')
    
    # Tracks generation
                                                  
    maskfilter = pe.Node(interface=niu.Function(input_names=['in_mask', 'filter_type', 'npass'], 
                                                  output_names=['out_file'], function=maskfilter), name='maskfilter')
    maskfilter.inputs.filter_type = 'erode'
    maskfilter.inputs.npass = 6
    
    dwi2response = pe.Node(interface=niu.Function(input_names=['in_dwi', 'in_bvecs', 'in_bvals', 'in_mask', 'lmax'], 
                                                  output_names=['response_file'], function=dwi2response), name='dwi2response')
    dwi2response.inputs.lmax = n_harmonics
    
    dwi2fod = pe.Node(interface=niu.Function(input_names=['in_dwi', 'in_bvecs', 'in_bvals', 'in_mask', 'lmax', 'response_file'], 
                                                  output_names=['SH_coeff_image'], function=dwi2fod), name='dwi2fod')
    dwi2fod.inputs.lmax = n_harmonics
    
    tckgen = pe.Node(interface=niu.Function(input_names=['algorithm', 'step_size', 'number_tracks', 'in_bvecs', 'in_bvals', 'white_matter_mask', 'SH_coeff_image'], 
                                            output_names=['out_tracks'], function=tckgen), name='tckgen')
    tckgen.inputs.algorithm = tractography_algorithm
    tckgen.inputs.step_size = step_size_tracking
    tckgen.inputs.number_tracks = nb_track
    
    # Outputnode and datasink
    
    outputnode = pe.Node(niu.IdentityInterface(fields=['tensor_file', 'FA_map', 'out_tracks', 'response_file','image_5tt','SH_coeff_image']), name='outputnode')
            
    datasink = pe.Node(nio.DataSink(), name='datasink')
    datasink.inputs.base_directory = op.join(datasink_directory, 'tractography/')
    
    # Workflow_construction    
    
    wf = pe.Workflow(name='Tractography')    

    wf.connect([(inputnode, dwi2tensor, [('in_dwi', 'in_dwi'),
                                         ('in_bvecs','in_bvecs'),
                                         ('in_bvals','in_bvals'),
                                         ('B0_mask','in_mask')])])
    wf.connect([(inputnode, tensor2fa_metric, [('B0_mask','in_mask')])])
    wf.connect([(dwi2tensor, tensor2fa_metric, [('tensor_file','tensor_file')])])
    wf.connect([(inputnode, maskfilter, [('B0_mask','in_mask')])])
    wf.connect([(inputnode, dwi2response, [('in_dwi', 'in_dwi'),
                                           ('in_bvecs','in_bvecs'),
                                           ('in_bvals','in_bvals')])])
    wf.connect([(maskfilter, dwi2response, [('out_file','in_mask')])])
    wf.connect([(inputnode, dwi2fod, [('in_dwi', 'in_dwi'),
                                      ('in_bvecs','in_bvecs'),
                                      ('in_bvals','in_bvals'),
                                      ('B0_mask','in_mask')])])                 
    wf.connect([(dwi2response, dwi2fod, [('response_file','response_file')])])
    wf.connect([(dwi2fod, tckgen, [('SH_coeff_image','SH_coeff_image')])])
    wf.connect([(inputnode, tckgen, [('in_bvecs','in_bvecs'),
                                      ('in_bvals','in_bvals')])])
    wf.connect([(inputnode, tckgen, [('white_matter_mask','white_matter_mask')])])
    wf.connect([(dwi2tensor, outputnode, [('tensor_file','tensor_file')])]) 
    wf.connect([(tensor2fa_metric, outputnode, [('FA_map','FA_map')])])
    wf.connect([(tckgen, outputnode, [('out_tracks','out_tracks')])])
    wf.connect([(dwi2response, outputnode, [('response_file','response_file')])])
    wf.connect([(dwi2fod, outputnode, [('SH_coeff_image','SH_coeff_image')])])
    wf.connect([(dwi2tensor, datasink, [('tensor_file','tensor_file')])]) 
    wf.connect([(tensor2fa_metric, datasink, [('FA_map','FA_map')])]) 
    wf.connect([(tckgen, datasink, [('out_tracks','out_tracks')])])
    wf.connect([(dwi2response, datasink, [('response_file','response_file')])])
    wf.connect([(dwi2fod, datasink, [('SH_coeff_image','SH_coeff_image')])])
    
    return wf
