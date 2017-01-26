import nipype.interfaces.spm.preprocess as spm
import nipype.interfaces.spm.utils as spmu
import nipype.pipeline.engine as pe
import nipype.interfaces.utility as niu
from nipype.interfaces.petpvc import PETPVC
from clinica.pipeline.pet.pet_utils import normalize_to_reference, create_mask
from clinica.pipeline.t1.t1_spm_utils import unzip_nii


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


    coreg_pet_t1 = pe.Node(spm.Coregister(), name='coreg_pet_t1')

    dartel_mni_reg = pe.Node(spm.DARTELNorm2MNI(), name='dartel_mni_reg')
    dartel_mni_reg.inputs.modulate = False
    dartel_mni_reg.inputs.fwhm = 0

    reslice = pe.Node(spmu.Reslice(), name='reslice')

    norm_to_ref = pe.Node(niu.Function(input_names=['pet_image', 'region_mask'],
                                       output_names=['suvr_pet_path'],
                                       function=normalize_to_reference),
                          name='norm_to_ref')

    outputnode = pe.Node(niu.IdentityInterface(fields=['suvr_pet_path']), name='outputnode')

    wf = pe.Workflow(name=name)
    if working_directory is not None:
        wf.base_dir = working_directory

    if pvc:
        inputnode = pe.Node(niu.IdentityInterface(
            fields=['pet_image', 't1_image_native', 'gm_image', 'wm_image', 'csf_image', 'fwhm_x', 'fwhm_y', 'fwhm_z',
                    'flow_fields', 'dartel_template', 'reference_mask']),
                            name='inputnode', mandatory_inputs=True)

        unzip_gm_image = pe.Node(niu.Function(input_names=['in_file'],
                                              output_names=['out_file'],
                                              function=unzip_nii),
                                 name='unzip_gm_image')

        unzip_wm_image = pe.Node(niu.Function(input_names=['in_file'],
                                              output_names=['out_file'],
                                              function=unzip_nii),
                                 name='unzip_wm_image')

        unzip_csf_image = pe.Node(niu.Function(input_names=['in_file'],
                                               output_names=['out_file'],
                                               function=unzip_nii),
                                  name='unzip_csf_image')

        mask = pe.Node(niu.Function(input_names=['c1', 'c2', 'c3'],
                                    output_names=['out_mask'],
                                    function=create_mask),
                       name='mask')

        petpvc = pe.Node(PETPVC(), name='pvc')
        petpvc.inputs.pvc = 'RBV'
        petpvc.inputs.out_file = 'pvc.nii'

        wf.connect([(inputnode, unzip_pet_image, [('pet_image', 'in_file')]),
                    (unzip_pet_image, coreg_pet_t1, [('out_file', 'source')]),
                    (inputnode, unzip_t1_image_native, [('t1_image_native', 'in_file')]),
                    (unzip_t1_image_native, coreg_pet_t1, [('out_file', 'target')]),

                    (inputnode, unzip_gm_image, [('gm_image', 'in_file')]),
                    (unzip_gm_image, mask, [('out_file', 'c1')]),
                    (inputnode, unzip_wm_image, [('wm_image', 'in_file')]),
                    (unzip_wm_image, mask, [('out_file', 'c2')]),
                    (inputnode, unzip_csf_image, [('csf_image', 'in_file')]),
                    (unzip_csf_image, mask, [('out_file', 'c3')]),
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
                    (mask, petpvc, [('out_mask', 'mask_file')]),
                    (petpvc, dartel_mni_reg, [('out_file', 'apply_to_files')]),
                    (dartel_mni_reg, reslice, [('normalized_files', 'space_defining')]),
                    (dartel_mni_reg, norm_to_ref, [('normalized_files', 'pet_image')]),
                    (reslice, norm_to_ref, [('out_file', 'region_mask')]),
                    (norm_to_ref, outputnode, [('suvr_pet_path', 'suvr_pet_path')])
                    ])
    else:
        inputnode = pe.Node(niu.IdentityInterface(
            fields=['pet_image', 't1_image_native', 'flow_fields', 'dartel_template', 'reference_mask']),
            name='inputnode', mandatory_inputs=True)

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
                    (norm_to_ref, outputnode, [('suvr_pet_path', 'suvr_pet_path')])
                    ])

    return wf
