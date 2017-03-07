import nipype.interfaces.spm.preprocess as spm
import nipype.interfaces.spm.utils as spmu
import nipype.pipeline.engine as pe
import nipype.interfaces.utility as niu
from nipype.interfaces.petpvc import PETPVC
from clinica.pipeline.pet.pet_utils import normalize_to_reference, create_binary_mask, apply_binary_mask, create_pvc_mask
from clinica.utils.io import zip_nii, unzip_nii


def pet_pipeline(working_directory=None,
                 name='pet_wf',
                 pvc=False):

    unzip_pet_image = pe.Node(niu.Function(input_names=['in_file'],
                                           output_names=['out_file'],
                                           function=unzip_nii),
                              name='unzip_pet_image')

    unzip_t1_image_native = pe.Node(niu.Function(input_names=['in_file'],
                                                 output_names=['out_file'],
                                                 function=unzip_nii),
                                    name='unzip_t1_image_native')

    unzip_flow_fields = pe.Node(niu.Function(input_names=['in_file'],
                                             output_names=['out_file'],
                                             function=unzip_nii),
                                name='unzip_flow_fields')

    unzip_dartel_template = pe.Node(niu.Function(input_names=['in_file'],
                                                 output_names=['out_file'],
                                                 function=unzip_nii),
                                    name='unzip_dartel_template')

    unzip_reference_mask = pe.Node(niu.Function(input_names=['in_file'],
                                                output_names=['out_file'],
                                                function=unzip_nii),
                                   name='unzip_reference_mask')

    unzip_mask_tissues = pe.MapNode(niu.Function(input_names=['in_file'],
                                                 output_names=['out_file'],
                                                 function=unzip_nii),
                                    name='unzip_mask_tissues',
                                    iterfield=['in_file'])

    coreg_pet_t1 = pe.Node(spm.Coregister(), name='coreg_pet_t1', copyfile=True)

    dartel_mni_reg = pe.Node(spm.DARTELNorm2MNI(), name='dartel_mni_reg', copyfile=True)
    dartel_mni_reg.inputs.modulate = False
    dartel_mni_reg.inputs.fwhm = 0

    reslice = pe.Node(spmu.Reslice(), name='reslice', copyfile=True)

    norm_to_ref = pe.Node(niu.Function(input_names=['pet_image', 'region_mask'],
                                       output_names=['suvr_pet_path'],
                                       function=normalize_to_reference),
                          name='norm_to_ref')

    binary_mask = pe.Node(niu.Function(input_names=['tissues'],
                                       output_names=['out_mask'],
                                       function=create_binary_mask),
                          name='binary_mask')

    apply_mask = pe.Node(niu.Function(input_names=['image', 'binary_mask'],
                                      output_names=['masked_image_path'],
                                      function=apply_binary_mask),
                         name='apply_mask')

    zip_pet_t1_native = pe.Node(niu.Function(input_names=['in_file'],
                                           output_names=['out_file'],
                                           function=zip_nii),
                              name='zip_pet_t1_native')

    zip_pet_mni = pe.Node(niu.Function(input_names=['in_file'],
                                                 output_names=['out_file'],
                                                 function=zip_nii),
                                    name='zip_pet_mni')

    zip_pet_suvr = pe.Node(niu.Function(input_names=['in_file'],
                                             output_names=['out_file'],
                                             function=zip_nii),
                                name='zip_pet_suvr')

    zip_binary_mask = pe.Node(niu.Function(input_names=['in_file'],
                                             output_names=['out_file'],
                                             function=zip_nii),
                                name='zip_binary_mask')

    zip_masked_pet_suvr = pe.Node(niu.Function(input_names=['in_file'],
                                               output_names=['out_file'],
                                               function=zip_nii),
                                  name='zip_masked_pet_suvr')

    wf = pe.Workflow(name=name)
    if working_directory is not None:
        wf.base_dir = working_directory

    if pvc:
        inputnode = pe.Node(niu.IdentityInterface(
            fields=['pet_image', 't1_image_native', 'mask_tissues', 'pvc_mask_tissues', 'fwhm_x', 'fwhm_y', 'fwhm_z',
                    'flow_fields', 'dartel_template', 'reference_mask']),
                            name='inputnode', mandatory_inputs=True)

        unzip_pvc_mask_tissues = pe.MapNode(niu.Function(input_names=['in_file'],
                                                     output_names=['out_file'],
                                                     function=unzip_nii),
                                        name='unzip_pvc_mask_tissues',
                                        iterfield=['in_file'])

        pvc_mask = pe.Node(niu.Function(input_names=['tissues'],
                                    output_names=['out_mask'],
                                    function=create_pvc_mask),
                       name='pvc_mask')

        petpvc = pe.Node(PETPVC(), name='pvc')
        petpvc.inputs.pvc = 'RBV'
        petpvc.inputs.out_file = 'pvc.nii'

        zip_pet_pvc = pe.Node(niu.Function(input_names=['in_file'],
                                           output_names=['out_file'],
                                           function=zip_nii),
                              name='zip_pet_pvc')

        outputnode = pe.Node(niu.IdentityInterface(fields=['pet_t1_native', 'pvc_pet', 'pvc_pet_mni',
                                                           'suvr_pvc_pet', 'binary_mask', 'masked_suvr_pvc_pet']),
                             name='outputnode')

        wf.connect([(inputnode, unzip_pet_image, [('pet_image', 'in_file')]),
                    (unzip_pet_image, coreg_pet_t1, [('out_file', 'source')]),
                    (inputnode, unzip_t1_image_native, [('t1_image_native', 'in_file')]),
                    (unzip_t1_image_native, coreg_pet_t1, [('out_file', 'target')]),

                    (inputnode, unzip_pvc_mask_tissues, [('pvc_mask_tissues', 'in_file')]),
                    (unzip_pvc_mask_tissues, pvc_mask, [('out_file', 'tissues')]),
                    (inputnode, petpvc, [('fwhm_x', 'fwhm_x'),
                                         ('fwhm_y', 'fwhm_y'),
                                         ('fwhm_z', 'fwhm_z')]),
                    (inputnode, unzip_flow_fields, [('flow_fields', 'in_file')]),
                    (unzip_flow_fields, dartel_mni_reg, [('out_file', 'flowfield_files')]),
                    (inputnode, unzip_dartel_template, [('dartel_template', 'in_file')]),
                    (unzip_dartel_template, dartel_mni_reg, [('out_file', 'template_file')]),
                    (inputnode, unzip_reference_mask, [('reference_mask', 'in_file')]),
                    (unzip_reference_mask, reslice, [('out_file', 'in_file')]),
                    (coreg_pet_t1, petpvc, [('coregistered_source', 'in_file')]),
                    (pvc_mask, petpvc, [('out_mask', 'mask_file')]),
                    (petpvc, dartel_mni_reg, [('out_file', 'apply_to_files')]),
                    (dartel_mni_reg, reslice, [('normalized_files', 'space_defining')]),
                    (dartel_mni_reg, norm_to_ref, [('normalized_files', 'pet_image')]),
                    (reslice, norm_to_ref, [('out_file', 'region_mask')]),

                    (coreg_pet_t1, zip_pet_t1_native, [('coregistered_source', 'in_file')]),
                    (petpvc, zip_pet_pvc, [('out_file', 'in_file')]),
                    (dartel_mni_reg, zip_pet_mni, [('normalized_files', 'in_file')]),
                    (norm_to_ref, zip_pet_suvr, [('suvr_pet_path', 'in_file')]),

                    (inputnode, unzip_mask_tissues, [('mask_tissues', 'in_file')]),
                    (unzip_mask_tissues, binary_mask, [('out_file', 'tissues')]),

                    (norm_to_ref, apply_mask, [('suvr_pet_path', 'image')]),
                    (binary_mask, apply_mask, [('out_mask', 'binary_mask')]),
                    (binary_mask, zip_binary_mask, [('out_mask', 'in_file')]),

                    (apply_mask, zip_masked_pet_suvr, [('masked_image_path', 'in_file')]),

                    (zip_pet_t1_native, outputnode, [('out_file', 'pet_t1_native')]),
                    (zip_pet_pvc, outputnode, [('out_file', 'pvc_pet')]),
                    (zip_pet_mni, outputnode, [('out_file', 'pvc_pet_mni')]),
                    (zip_pet_suvr, outputnode, [('out_file', 'suvr_pvc_pet')]),
                    (zip_binary_mask, outputnode, [('out_file', 'binary_mask')]),
                    (zip_masked_pet_suvr, outputnode, [('out_file', 'masked_suvr_pvc_pet')])
                    ])
    else:
        inputnode = pe.Node(niu.IdentityInterface(
            fields=['pet_image', 't1_image_native', 'mask_tissues', 'flow_fields', 'dartel_template', 'reference_mask']),
            name='inputnode', mandatory_inputs=True)

        outputnode = pe.Node(niu.IdentityInterface(fields=['pet_t1_native', 'pet_mni', 'suvr_pet', 'binary_mask', 'masked_suvr_pet']),
                             name='outputnode')

        wf.connect([(inputnode, unzip_pet_image, [('pet_image', 'in_file')]),
                    (unzip_pet_image, coreg_pet_t1, [('out_file', 'source')]),
                    (inputnode, unzip_t1_image_native, [('t1_image_native', 'in_file')]),
                    (unzip_t1_image_native, coreg_pet_t1, [('out_file', 'target')]),
                    (inputnode, unzip_flow_fields, [('flow_fields', 'in_file')]),
                    (unzip_flow_fields, dartel_mni_reg, [('out_file', 'flowfield_files')]),
                    (inputnode, unzip_dartel_template, [('dartel_template', 'in_file')]),
                    (unzip_dartel_template, dartel_mni_reg, [('out_file', 'template_file')]),
                    (inputnode, unzip_reference_mask, [('reference_mask', 'in_file')]),
                    (unzip_reference_mask, reslice, [('out_file', 'in_file')]),
                    (coreg_pet_t1, dartel_mni_reg, [('coregistered_source', 'apply_to_files')]),
                    (dartel_mni_reg, reslice, [('normalized_files', 'space_defining')]),
                    (dartel_mni_reg, norm_to_ref, [('normalized_files', 'pet_image')]),
                    (reslice, norm_to_ref, [('out_file', 'region_mask')]),

                    (coreg_pet_t1, zip_pet_t1_native, [('coregistered_source', 'in_file')]),
                    (dartel_mni_reg, zip_pet_mni, [('normalized_files', 'in_file')]),
                    (norm_to_ref, zip_pet_suvr, [('suvr_pet_path', 'in_file')]),

                    (inputnode, unzip_mask_tissues, [('mask_tissues', 'in_file')]),
                    (unzip_mask_tissues, binary_mask, [('out_file', 'tissues')]),

                    (norm_to_ref, apply_mask, [('suvr_pet_path', 'image')]),
                    (binary_mask, apply_mask, [('out_mask', 'binary_mask')]),
                    (binary_mask, zip_binary_mask, [('out_mask', 'in_file')]),
                    (apply_mask, zip_masked_pet_suvr, [('masked_image_path', 'in_file')]),

                    (zip_pet_t1_native, outputnode, [('out_file', 'pet_t1_native')]),
                    (zip_pet_mni, outputnode, [('out_file', 'pet_mni')]),
                    (zip_pet_suvr, outputnode, [('out_file', 'suvr_pet')]),
                    (zip_binary_mask, outputnode, [('out_file', 'binary_mask')]),
                    (zip_masked_pet_suvr, outputnode, [('out_file', 'masked_suvr_pet')])
                    ])

    return wf
