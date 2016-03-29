# -*- coding: utf-8 -*-
"""
Created on Wed Mar 16 16:23:47 2016

@author: jacquemont
"""

def launch(in_dwi, in_T1, in_bvals, in_bvecs, working_directory, datasink_directory):
    """
    Create and run a high level pipeline to preprocess the DWI Images :
        - Preparation of the dataset
        - Correction for Head Motion 
        - Correction for Eddy Currents 
        - Correction for EPI susceptibility induced distortions without the field map
        - Bias field correction
    The outputs presented are tipically outputs necessary for further tractography.
    
    Inputs
    ---------
    in_dwi : STRING
      Path to the DWI image.
    in_T1: STRING
      Path to the T1 image.
    in_bvals: STRING
      Path to the b-vals text file.
    in_bvecs: STRING
      Path to the b-vecs text file.
    working_directory : STRING
      Directory to use as tmp for all the temporary files generated by the workflow.
    datasink_directory : STRING
      Base directory of the datasink.
    
    Outputs
    ----------
        DWI_hmc_ecc_epi_bias_corrected - DWI corrected for Head motion, Eddy currents, EPI susceptibility induced distortions and bias field
        out_bvecs - updated and corrected gradient vectors table
        out_bvals - updated gradient values table
        mask_b0 - Binary mask obtained from the average of the B0 images    
    
    """
    
    import nipype.interfaces.io as nio
    import nipype.interfaces.utility as niu
    import nipype.pipeline.engine as pe
    import os.path as op
    import clinica.pipeline.preprocessing.DWI_corrections as predifcorrect

# Inputs existence checking

    inputs=[in_dwi, in_T1, in_bvals, in_bvecs, working_directory, datasink_directory]     
        
    for input_file in inputs:
        if not op.exists(input_file):
            raise IOError('file {} does not exist'.format(input_file))
    
    datasource = pe.Node(interface=nio.DataGrabber(infields=[], outfields=['dwi_image','bvectors_directions','bvalues','T1_image']), name='datasource')
    datasource.inputs.template = '*'
    datasource.inputs.field_template = dict(dwi_image= in_dwi,
                                            bvalues=in_bvals,
                                            bvectors_directions= in_bvecs,
                                            T1_image= in_T1)
    datasource.inputs.template_args = dict(dwi_image=[[]],
                                           bvalues=[[]],
                                           bvectors_directions=[[]],
                                           T1_image=[[]])
    datasource.inputs.sort_filelist = True
    
    inputnode = pe.Node(interface=niu.IdentityInterface(fields=["dwi_image", "bvectors_directions", "bvalues", 'T1_image']), name="inputnode")
    
    pre = predifcorrect.prepare_data(datasink_directory)
    
    hmc = predifcorrect.hmc_pipeline(datasink_directory)
    
    ecc = predifcorrect.ecc_pipeline(datasink_directory)

    epi = predifcorrect.epi_pipeline(datasink_directory)

    bias = predifcorrect.remove_bias(datasink_directory)
    
    aac = predifcorrect.apply_all_corrections(datasink_directory)
    
    datasink = pe.Node(nio.DataSink(), name='datasink_tracto')
    datasink.inputs.base_directory = op.join(datasink_directory, 'Outputs_for_Tractography/')
    
    wf = pe.Workflow(name='preprocess')
    wf.base_dir = working_directory
    
    wf.connect([(datasource, inputnode, [('dwi_image','dwi_image'), ('bvalues','bvalues'), ('bvectors_directions','bvectors_directions'), ('T1_image','T1_image')])])
    wf.connect([(inputnode, pre, [('dwi_image', 'inputnode.dwi_image'),
                                  ('bvalues', 'inputnode.bvalues'),
                                  ('bvectors_directions', 'inputnode.bvectors_directions')])])
    wf.connect([(pre, hmc,[('outputnode.dwi_b0_merge','inputnode.in_file'), ('outputnode.out_bvals','inputnode.in_bval'), ('outputnode.out_bvecs','inputnode.in_bvec')])])
    wf.connect([(pre, hmc, [('outputnode.mask_b0','inputnode.in_mask')])])
    wf.connect([(hmc, ecc, [('outputnode.out_xfms','inputnode.in_xfms'),('outputnode.out_file','inputnode.in_file')])])
    wf.connect([(pre, ecc, [('outputnode.out_bvals','inputnode.in_bval')])])
    wf.connect([(pre, ecc, [('outputnode.mask_b0','inputnode.in_mask')])])
    wf.connect([(ecc, epi, [('outputnode.out_file','inputnode.DWI')])])
    wf.connect([(inputnode, epi, [('T1_image','inputnode.T1')])])
    wf.connect([(hmc, epi, [('outputnode.out_bvec','inputnode.bvec')])])
    wf.connect([(pre, aac, [('outputnode.dwi_b0_merge', 'inputnode.in_dwi')])])
    wf.connect([(hmc, aac, [('outputnode.out_xfms', 'inputnode.in_hmc')])])
    wf.connect([(ecc, aac, [('outputnode.out_xfms', 'inputnode.in_ecc')])])
    wf.connect([(epi, aac, [('outputnode.out_warp', 'inputnode.in_epi')])])
    wf.connect([(inputnode, aac, [('T1_image','inputnode.T1')])])
    
    wf.connect([(aac, bias, [('outputnode.out_file','inputnode.in_file')])])
    
    wf.connect([(bias, datasink, [('outputnode.out_file','DWI_hmc_ecc_epi_bias_corrected')])])
    wf.connect([(epi, datasink, [('outputnode.out_bvec', 'out_bvecs')])])
    wf.connect([(pre, datasink, [('outputnode.out_bvals','out_bvals')])])
    wf.connect([(bias, datasink, [('outputnode.b0_mask','b0_mask')])])
    
    wf.run()
