"""PET Preprocess Volume - Clinica Pipeline.
This file has been generated automatically by the `clinica generate template`
command line tool. See here for more details: https://gitlab.icm-institute.org/aramislab/clinica/wikis/docs/InteractingWithClinica.
"""

import clinica.pipeline.engine as cpe


class PETPreprocessVolume(cpe.Pipeline):
    """PET Preprocess Volume SHORT DESCRIPTION.

    Warnings:
        - A WARNING.

    Todos:
        - [x] A FILLED TODO ITEM.
        - [ ] AN ON-GOING TODO ITEM.

    Args:
        input_dir: A BIDS directory.
        output_dir: An empty output directory where CAPS structured data will be written.
        subjects_sessions_list: The Subjects-Sessions list file (in .tsv format).

    Returns:
        A clinica pipeline object containing the PET Preprocess Volume pipeline.

    Raises:


    Example:
        >>> from pet_preprocess_volume import PETPreprocessVolume
        >>> pipeline = PETPreprocessVolume('myGroup', '~/MYDATASET_BIDS', '~/MYDATASET_CAPS')
        >>> pipeline.parameters = {
        >>>     # ...
        >>> }
        >>> pipeline.base_dir = '/tmp/'
        >>> pipeline.run()
    """

    def __init__(self, bids_directory=None, caps_directory=None, tsv_file=None, name=None, group_id='default', fwhm_tsv=None):
        super(PETPreprocessVolume, self).__init__(bids_directory, caps_directory, tsv_file, name)

        if not group_id.isalnum():
            raise ValueError('Not valid group_id value. It must be composed only by letters and/or numbers')

        from pandas.io.parsers import read_csv

        self._group_id = group_id
        self._suvr_region = ''
        self._fwhm = None
        self._apply_pvc = False

        if fwhm_tsv is not None:
            fwhm_df = read_csv(fwhm_tsv, sep='\t')

            if fwhm_df.shape[0] != len(self.subjects):
                raise ValueError('The number of rows in fwhm_tsv file must match the number of subject-session pairs.')

            fwhm_x = list(fwhm_df.fwhm_x)
            fwhm_y = list(fwhm_df.fwhm_y)
            fwhm_z = list(fwhm_df.fwhm_z)
            self._fwhm = [[fwhm_x[i], fwhm_y[i], fwhm_z[i]] for i in range(len(self.subjects))]
            self._apply_pvc = True

        # Default parameters
        self._parameters = {'pet_type': 'fdg',
                            'mask_tissues': [1, 2, 3],
                            'mask_tissues_modulation': 'on',
                            'mask_threshold': 0.3,
                            'pvc_mask_tissues': [1, 2, 3],
                            'atlas_list': ['AAL2', 'LPBA40', 'Neuromorphometrics', 'AICHA', 'Hammers']
                            }

    def check_custom_dependencies(self):
        """Check dependencies that can not be listed in the `info.json` file.
        """
        pass

    def get_input_fields(self):
        """Specify the list of possible inputs of this pipeline.

        Returns:
            A list of (string) input fields name.
        """

        return ['pet_image',
                't1_image_native',
                'mask_tissues',
                'pvc_mask_tissues',
                'fwhm',
                'flow_fields',
                'dartel_template',
                'reference_mask']

    def get_output_fields(self):
        """Specify the list of possible outputs of this pipeline.

        Returns:
            A list of (string) output fields name.
        """

        return ['pet_t1_native',
                'pet_mni',
                'pet_suvr',
                'binary_mask',
                'pet_suvr_masked',
                'pet_pvc',
                'pet_pvc_mni',
                'pet_pvc_suvr',
                'pet_pvc_suvr_masked',
                'atlas_statistics',
                'pvc_atlas_statistics'
                ]

    def build_input_node(self):
        """Build and connect an input node to the pipeline.
        """

        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        import nipype.interfaces.io as nio
        import clinica.pipeline.pet_preprocess_volume.pet_preprocess_volume_utils as utils
        import os
        from os.path import join

        iterables_fwhm = self._fwhm
        if not self._apply_pvc:
            iterables_fwhm = [[]] * len(self.subjects)

        iterables_node = npe.Node(name="LoadingCLIArguments",
                                  interface=nutil.IdentityInterface(
                                      fields=['subject_id', 'session_id', 'fwhm'],
                                      mandatory_inputs=True))
        iterables_node.iterables = [('subject_id', self.subjects),
                                    ('session_id', self.sessions),
                                    ('fwhm', iterables_fwhm)]
        iterables_node.synchronize = True

        # PET DataGrabber
        # ===============
        pet_bids_reader = npe.Node(nio.DataGrabber(infields=['subject_id', 'session', 'subject_repeat', 'session_repeat'],
                                                   outfields=['out_files']), name='pet_bids_reader')
        pet_bids_reader.inputs.base_directory = self.bids_directory
        pet_bids_reader.inputs.template = '%s/%s/pet/%s_%s_task-rest_acq-' + self.parameters['pet_type'] + '_pet.nii*'
        pet_bids_reader.inputs.sort_filelist = False

        # Native T1 DataGrabber
        # ======================
        t1_bids_reader = npe.Node(nio.DataGrabber(infields=['subject_id', 'session', 'subject_repeat', 'session_repeat'],
                                                  outfields=['out_files']), name='t1_bids_reader')
        t1_bids_reader.inputs.base_directory = self.bids_directory
        t1_bids_reader.inputs.template = '%s/%s/anat/%s_%s_T1w.nii*'
        t1_bids_reader.inputs.sort_filelist = False

        # Flow Fields DataGrabber
        # ========================
        flowfields_caps_reader = npe.Node(nio.DataGrabber(infields=['subject_id', 'session', 'subject_repeat', 'session_repeat'],
                                                   outfields=['out_files']), name='flowfields_caps_reader')
        flowfields_caps_reader.inputs.base_directory = join(self.caps_directory, 'subjects')
        flowfields_caps_reader.inputs.template = '%s/%s/t1/spm/dartel/group-' + self._group_id + '/%s_%s_T1w_target-' + self._group_id + '_transformation-forward_deformation.nii*'
        flowfields_caps_reader.inputs.sort_filelist = False


        # Dartel Template DataGrabber
        # ============================
        template_caps_reader = npe.Node(nio.DataGrabber(outfields=['out_files']),
                                        name="template_caps_reader")
        template_caps_reader.inputs.base_directory = self.caps_directory
        template_caps_reader.inputs.template = 'groups/group-' + self._group_id + '/t1/group-' + self._group_id + '_template.nii*'
        template_caps_reader.inputs.sort_filelist = False

        # Reference Mask DataGrabber
        # ===========================
        clinica_home = os.environ.get("CLINICA_HOME")
        reference_mask = npe.Node(nio.DataGrabber(outfields=['out_files']), name='reference_mask')
        reference_mask.inputs.base_directory = join(clinica_home, 'clinica', 'resources', 'masks')
        reference_mask.inputs.sort_filelist = False
        # TODO DIFFERENT PET TYPES TO PROCESS
        if self.parameters['pet_type'] == 'fdg' or self.parameters['pet_type'] == 'FDG':
            reference_mask.inputs.template = 'region-pons_eroded-6mm_mask.nii*'
            self._suvr_region = 'pons'
        elif self.parameters['pet_type'] == 'av45':
            reference_mask.inputs.template = 'region-cerebellumPons_eroded-6mm_mask.nii*'
            self._suvr_region = 'cerebellumPons'
        else:
            raise NotImplementedError('Unknown type of PET image. We currently accept as input only "fdg" or "av45" as values.')

        # Tissues DataGrabber
        # ====================
        tissue_names = {1: 'graymatter',
                        2: 'whitematter',
                        3: 'csf',
                        4: 'bone',
                        5: 'softtissue',
                        6: 'background'
                        }

        tissues_caps_reader = npe.Node(nio.DataGrabber(infields=['subject_id', 'session', 'subject_repeat', 'session_repeat', 'tissues'],
                                                       outfields=['out_files']), name='tissues_caps_reader')
        tissues_caps_reader.inputs.base_directory = join(self.caps_directory, 'subjects')
        tissues_caps_reader.inputs.template = '%s/%s/t1/spm/dartel/group-' + self._group_id + '/%s_%s_T1w_segm-%s_space-Ixi549Space_modulated-' + self.parameters['mask_tissues_modulation'] + '_probability.nii*'
        tissues_caps_reader.inputs.tissues = [tissue_names[t] for t in self.parameters['mask_tissues']]
        tissues_caps_reader.inputs.sort_filelist = False

        n_tissues = len(self.parameters['mask_tissues'])

        self.connect([(iterables_node, pet_bids_reader, [('subject_id', 'subject_id'),
                                                         ('session_id', 'session'),
                                                         ('subject_id', 'subject_repeat'),
                                                         ('session_id', 'session_repeat')]),
                      (iterables_node, t1_bids_reader, [('subject_id', 'subject_id'),
                                                        ('session_id', 'session'),
                                                        ('subject_id', 'subject_repeat'),
                                                        ('session_id', 'session_repeat')]),
                      (iterables_node, tissues_caps_reader, [(('subject_id', utils.expand_into_list, n_tissues), 'subject_id'),
                                                             (('session_id', utils.expand_into_list, n_tissues), 'session'),
                                                             (('subject_id', utils.expand_into_list, n_tissues), 'subject_repeat'),
                                                             (('session_id', utils.expand_into_list, n_tissues), 'session_repeat')]),
                      (iterables_node, flowfields_caps_reader, [('subject_id', 'subject_id'),
                                                                ('session_id', 'session'),
                                                                ('subject_id', 'subject_repeat'),
                                                                ('session_id', 'session_repeat')]),
                      (pet_bids_reader, self.input_node, [('out_files', 'pet_image')]),
                      (t1_bids_reader, self.input_node, [('out_files', 't1_image_native')]),
                      (tissues_caps_reader, self.input_node, [('out_files', 'mask_tissues')]),
                      (flowfields_caps_reader, self.input_node, [('out_files', 'flow_fields')]),
                      (template_caps_reader, self.input_node, [('out_files', 'dartel_template')]),
                      (reference_mask, self.input_node, [('out_files', 'reference_mask')]),
                      (iterables_node, self.input_node, [('fwhm', 'fwhm')])
                      ])
        if self._apply_pvc:
            pvc_tissues_caps_reader = npe.Node(
                nio.DataGrabber(infields=['subject_id', 'session', 'subject_repeat', 'session_repeat', 'tissues'],
                                outfields=['out_files']), name='pvc_tissues_caps_reader')
            pvc_tissues_caps_reader.inputs.base_directory = join(self.caps_directory, 'subjects')
            pvc_tissues_caps_reader.inputs.template = '%s/%s/t1/spm/segmentation/native_space/%s_%s_T1w_segm-%s_probability.nii*'
            pvc_tissues_caps_reader.inputs.tissues = [tissue_names[t] for t in self.parameters['pvc_mask_tissues']]
            pvc_tissues_caps_reader.inputs.sort_filelist = False

            n_pvc_tissues = len(self.parameters['pvc_mask_tissues'])

            self.connect([(iterables_node, pvc_tissues_caps_reader, [(('subject_id', utils.expand_into_list, n_pvc_tissues), 'subject_id'),
                                                                     (('session_id', utils.expand_into_list, n_pvc_tissues), 'session'),
                                                                     (('subject_id', utils.expand_into_list, n_pvc_tissues), 'subject_repeat'),
                                                                     (('session_id', utils.expand_into_list, n_pvc_tissues), 'session_repeat')]),
                          (pvc_tissues_caps_reader, self.input_node, [('out_files', 'pvc_mask_tissues')])
                          ])
        else:
            self.input_node.inputs.pvc_mask_tissues = []

    def build_output_node(self):
        """Build and connect an output node to the pipeline.
        """

        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        import nipype.interfaces.io as nio
        from clinica.utils.io import zip_nii
        from os.path import join
        import re
        import clinica.pipeline.pet_preprocess_volume.pet_preprocess_volume_utils as utils

        # Find container path from pet filename
        # =====================================
        container_path = npe.Node(nutil.Function(input_names=['pet_filename'],
                                                 output_names=['container'],
                                                 function=utils.pet_container_from_filename),
                                  name='container_path')
        container_path.inputs.threshold = self.parameters['mask_threshold']

        # Writing all images into CAPS
        # ============================
        write_images_node = npe.Node(name='write_caps_node',
                                     interface=nio.DataSink())
        write_images_node.inputs.base_directory = self.caps_directory
        write_images_node.inputs.parameterization = False
        write_images_node.inputs.regexp_substitutions = [
            (r'(.*/)pet_t1_native/r(sub-.*)(\.nii(\.gz)?)$', r'\1\2_space-T1w_pet\3'),
            (r'(.*/)pet_pvc/pvc-rbv_r(sub-.*)(\.nii(\.gz)?)$', r'\1\2_space-T1w_pvc-rbv_pet\3'),
            (r'(.*/)pet_mni/wr(sub-.*)(\.nii(\.gz)?)$', r'\1\2_space-Ixi549Space_pet\3'),
            (r'(.*/)pet_pvc_mni/wpvc-rbv_r(sub-.*)(\.nii(\.gz)?)$', r'\1\2_space-Ixi549Space_pvc-rbv_pet\3'),
            (r'(.*/)pet_suvr/suvr_wr(sub-.*)(\.nii(\.gz)?)$', r'\1\2_space-Ixi549Space_suvr-' + re.escape(self._suvr_region) + r'_pet\3'),
            (r'(.*/)pet_pvc_suvr/suvr_wpvc-rbv_r(sub-.*)(\.nii(\.gz)?)$', r'\1\2_space-Ixi549Space_pvc-rbv_suvr-' + re.escape(self._suvr_region) + r'_pet\3'),
            (r'(.*/)pet_suvr_masked/masked_suvr_wr(sub-.*)(\.nii(\.gz)?)$', r'\1\2_space-Ixi549Space_suvr-' + re.escape(self._suvr_region) + r'_mask-brain_pet\3'),
            (r'(.*/)pet_pvc_suvr_masked/masked_suvr_wpvc-rbv_r(sub-.*)(\.nii(\.gz)?)$', r'\1\2_space-Ixi549Space_pvc-rbv_suvr-' + re.escape(self._suvr_region) + r'_mask-brain_pet\3'),
            (r'(.*/)binary_mask/(sub-.*_T1w_).*(space-[a-zA-Z0-9]+).*(_brainmask\.nii(\.gz)?)$', r'\1\2\3\4')
        ]

        # Writing atlas statistics into CAPS
        # ==================================
        write_atlas_node = npe.Node(name='write_atlas_node',
                                    interface=nio.DataSink())
        write_atlas_node.inputs.base_directory = self.caps_directory
        write_atlas_node.inputs.parameterization = False
        write_atlas_node.inputs.regexp_substitutions = [
            (r'(.*/atlas_statistics/)suvr_wr(sub-.*)(_statistics\.tsv)$', r'\1\2' + r'_suvr-' +
             re.escape(self._suvr_region) + r'\3'),
            (r'(.*/)pvc_(atlas_statistics/)suvr_wpvc-rbv_r(sub-.*)(_statistics\.tsv)$', r'\1\2\3' + r'_pvc-rbv_suvr-' +
             re.escape(self._suvr_region) + r'\4')
        ]

        self.connect([(self.input_node, container_path, [('pet_image', 'pet_filename')]),
                      (container_path, write_images_node, [(('container', join, 'group-' + self._group_id), 'container')]),
                      (self.output_node, write_images_node, [(('pet_t1_native', zip_nii, True), 'pet_t1_native'),
                                                             (('pet_mni', zip_nii, True), 'pet_mni'),
                                                             (('pet_suvr', zip_nii, True), 'pet_suvr'),
                                                             (('binary_mask', zip_nii, True), 'binary_mask'),
                                                             (('pet_suvr_masked', zip_nii, True), 'pet_suvr_masked'),
                                                             (('pet_pvc', zip_nii, True), 'pet_pvc'),
                                                             (('pet_pvc_mni', zip_nii, True), 'pet_pvc_mni'),
                                                             (('pet_pvc_suvr', zip_nii, True), 'pet_pvc_suvr'),
                                                             (('pet_pvc_suvr_masked', zip_nii, True), 'pet_pvc_suvr_masked')]),
                      (container_path, write_atlas_node, [('container', 'container')]),
                      (self.output_node, write_atlas_node, [('atlas_statistics', 'atlas_statistics'),
                                                            ('pvc_atlas_statistics', 'pvc_atlas_statistics')])
                      ])

    def build_core_nodes(self):
        """Build and connect the core nodes of the pipeline.
        """

        import nipype.interfaces.spm as spm
        import nipype.interfaces.spm.utils as spmutils
        from nipype.interfaces.petpvc import PETPVC
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        from clinica.utils.io import unzip_nii
        import clinica.pipeline.pet_preprocess_volume.pet_preprocess_volume_utils as utils

        # Unzipping
        # =========
        unzip_pet_image = npe.Node(nutil.Function(input_names=['in_file'],
                                                  output_names=['out_file'],
                                                  function=unzip_nii),
                                   name='unzip_pet_image')

        unzip_t1_image_native = npe.Node(nutil.Function(input_names=['in_file'],
                                                        output_names=['out_file'],
                                                        function=unzip_nii),
                                         name='unzip_t1_image_native')

        unzip_flow_fields = npe.Node(nutil.Function(input_names=['in_file'],
                                                    output_names=['out_file'],
                                                    function=unzip_nii),
                                     name='unzip_flow_fields')

        unzip_dartel_template = npe.Node(nutil.Function(input_names=['in_file'],
                                                        output_names=['out_file'],
                                                        function=unzip_nii),
                                         name='unzip_dartel_template')

        unzip_reference_mask = npe.Node(nutil.Function(input_names=['in_file'],
                                                       output_names=['out_file'],
                                                       function=unzip_nii),
                                        name='unzip_reference_mask')

        unzip_mask_tissues = npe.MapNode(nutil.Function(input_names=['in_file'],
                                                        output_names=['out_file'],
                                                        function=unzip_nii),
                                         name='unzip_mask_tissues',
                                         iterfield=['in_file'])

        # Coregister PET into T1 native space
        # ===================================
        coreg_pet_t1 = npe.Node(spm.Coregister(), name='coreg_pet_t1')

        # Spatially normalize PET into MNI
        # ================================
        dartel_mni_reg = npe.Node(spm.DARTELNorm2MNI(), name='dartel_mni_reg')
        dartel_mni_reg.inputs.modulate = False
        dartel_mni_reg.inputs.fwhm = 0

        # Reslice reference region mask into PET
        # ======================================
        reslice = npe.Node(spmutils.Reslice(), name='reslice')

        # Normalize PET values according to reference region
        # ==================================================
        norm_to_ref = npe.Node(nutil.Function(input_names=['pet_image', 'region_mask'],
                                              output_names=['suvr_pet_path'],
                                              function=utils.normalize_to_reference),
                               name='norm_to_ref')

        # Create binary mask from segmented tissues
        # =========================================
        binary_mask = npe.Node(nutil.Function(input_names=['tissues', 'threshold'],
                                              output_names=['out_mask'],
                                              function=utils.create_binary_mask),
                               name='binary_mask')
        binary_mask.inputs.threshold = self.parameters['mask_threshold']

        # Mask PET image
        # ==============
        apply_mask = npe.Node(nutil.Function(input_names=['image', 'binary_mask'],
                                             output_names=['masked_image_path'],
                                             function=utils.apply_binary_mask),
                              name='apply_mask')

        # Atlas Statistics
        # ================
        atlas_stats_node = npe.MapNode(nutil.Function(input_names=['in_image',
                                                                   'in_atlas_list'],
                                                      output_names=['atlas_statistics'],
                                                      function=utils.atlas_statistics),
                                       name='atlas_stats_node',
                                       iterfield=['in_image'])
        atlas_stats_node.inputs.in_atlas_list = self.parameters['atlas_list']

        # Connection
        # ==========
        self.connect([(self.input_node, unzip_pet_image, [('pet_image', 'in_file')]),
                      (self.input_node, unzip_t1_image_native, [('t1_image_native', 'in_file')]),
                      (self.input_node, unzip_flow_fields, [('flow_fields', 'in_file')]),
                      (self.input_node, unzip_dartel_template, [('dartel_template', 'in_file')]),
                      (self.input_node, unzip_reference_mask, [('reference_mask', 'in_file')]),
                      (self.input_node, unzip_mask_tissues, [('mask_tissues', 'in_file')]),

                      (unzip_pet_image, coreg_pet_t1, [('out_file', 'source')]),
                      (unzip_t1_image_native, coreg_pet_t1, [('out_file', 'target')]),
                      (unzip_flow_fields, dartel_mni_reg, [('out_file', 'flowfield_files')]),
                      (unzip_dartel_template, dartel_mni_reg, [('out_file', 'template_file')]),
                      (unzip_reference_mask, reslice, [('out_file', 'in_file')]),
                      (unzip_mask_tissues, binary_mask, [('out_file', 'tissues')]),

                      (coreg_pet_t1, dartel_mni_reg, [('coregistered_source', 'apply_to_files')]),
                      (dartel_mni_reg, reslice, [('normalized_files', 'space_defining')]),
                      (dartel_mni_reg, norm_to_ref, [('normalized_files', 'pet_image')]),
                      (reslice, norm_to_ref, [('out_file', 'region_mask')]),
                      (norm_to_ref, apply_mask, [('suvr_pet_path', 'image')]),
                      (binary_mask, apply_mask, [('out_mask', 'binary_mask')]),
                      (norm_to_ref, atlas_stats_node, [('suvr_pet_path', 'in_image')]),

                      (coreg_pet_t1, self.output_node, [('coregistered_source', 'pet_t1_native')]),
                      (dartel_mni_reg, self.output_node, [('normalized_files', 'pet_mni')]),
                      (norm_to_ref, self.output_node, [('suvr_pet_path', 'pet_suvr')]),
                      (binary_mask, self.output_node, [('out_mask', 'binary_mask')]),
                      (apply_mask, self.output_node, [('masked_image_path', 'pet_suvr_masked')]),
                      (atlas_stats_node, self.output_node, [('atlas_statistics', 'atlas_statistics')])
                      ])

        # PVC
        # ==========
        if self._apply_pvc:
            # Unzipping
            # =========
            unzip_pvc_mask_tissues = npe.MapNode(nutil.Function(input_names=['in_file'],
                                                                output_names=['out_file'],
                                                                function=unzip_nii),
                                                 name='unzip_pvc_mask_tissues',
                                                 iterfield=['in_file'])

            # Creating Mask to use in PVC
            # ===========================
            pvc_mask = npe.Node(nutil.Function(input_names=['tissues'],
                                               output_names=['out_mask'],
                                               function=utils.create_pvc_mask),
                                name='pvc_mask')
            # PET PVC
            # =======
            petpvc = npe.Node(PETPVC(), name='pvc')
            petpvc.inputs.pvc = 'RBV'
            petpvc.inputs.out_file = 'pvc.nii'

            # Spatially normalize PET into MNI
            # ================================
            dartel_mni_reg_pvc = npe.Node(spm.DARTELNorm2MNI(), name='dartel_mni_reg_pvc')
            dartel_mni_reg_pvc.inputs.modulate = False
            dartel_mni_reg_pvc.inputs.fwhm = 0

            # Reslice reference region mask into PET
            # ======================================
            reslice_pvc = npe.Node(spmutils.Reslice(), name='reslice_pvc')

            # Normalize PET values according to reference region
            # ==================================================
            norm_to_ref_pvc = npe.Node(nutil.Function(input_names=['pet_image', 'region_mask'],
                                                  output_names=['suvr_pet_path'],
                                                  function=utils.normalize_to_reference),
                                   name='norm_to_ref_pvc')

            # Mask PET image
            # ==============
            apply_mask_pvc = npe.Node(nutil.Function(input_names=['image', 'binary_mask'],
                                                     output_names=['masked_image_path'],
                                                     function=utils.apply_binary_mask),
                                      name='apply_mask_pvc')
            # Atlas Statistics
            # ================
            atlas_stats_pvc = npe.MapNode(nutil.Function(input_names=['in_image',
                                                                      'in_atlas_list'],
                                                         output_names=['atlas_statistics'],
                                                         function=utils.atlas_statistics),
                                          name='atlas_stats_pvc',
                                          iterfield=['in_image'])
            atlas_stats_pvc.inputs.in_atlas_list = self.parameters['atlas_list']

            # Connection
            # ==========
            self.connect([(self.input_node, unzip_pvc_mask_tissues, [('pvc_mask_tissues', 'in_file')]),
                          (unzip_pvc_mask_tissues, pvc_mask, [('out_file', 'tissues')]),
                          (unzip_flow_fields, dartel_mni_reg_pvc, [('out_file', 'flowfield_files')]),
                          (unzip_dartel_template, dartel_mni_reg_pvc, [('out_file', 'template_file')]),
                          (unzip_reference_mask, reslice_pvc, [('out_file', 'in_file')]),
                          (coreg_pet_t1, petpvc, [('coregistered_source', 'in_file'),
                                                  (('coregistered_source', utils.pet_pvc_name, 'RBV'), 'out_file')]),
                          (pvc_mask, petpvc, [('out_mask', 'mask_file')]),
                          (self.input_node, petpvc, [(('fwhm', utils.get_from_list, 0), 'fwhm_x'),
                                                     (('fwhm', utils.get_from_list, 1), 'fwhm_y'),
                                                     (('fwhm', utils.get_from_list, 2), 'fwhm_z')]),
                          (petpvc, dartel_mni_reg_pvc, [('out_file', 'apply_to_files')]),
                          (dartel_mni_reg_pvc, reslice_pvc, [('normalized_files', 'space_defining')]),
                          (dartel_mni_reg_pvc, norm_to_ref_pvc, [('normalized_files', 'pet_image')]),
                          (reslice_pvc, norm_to_ref_pvc, [('out_file', 'region_mask')]),

                          (norm_to_ref_pvc, apply_mask_pvc, [('suvr_pet_path', 'image')]),
                          (binary_mask, apply_mask_pvc, [('out_mask', 'binary_mask')]),
                          (norm_to_ref_pvc, atlas_stats_pvc, [('suvr_pet_path', 'in_image')]),

                          (petpvc, self.output_node, [('out_file', 'pet_pvc')]),
                          (dartel_mni_reg_pvc, self.output_node, [('normalized_files', 'pet_pvc_mni')]),
                          (norm_to_ref_pvc, self.output_node, [('suvr_pet_path', 'pet_pvc_suvr')]),
                          (apply_mask_pvc, self.output_node, [('masked_image_path', 'pet_pvc_suvr_masked')]),
                          (atlas_stats_pvc, self.output_node, [('atlas_statistics', 'pvc_atlas_statistics')])
                          ])
        else:
            self.output_node.inputs.pet_pvc = []
            self.output_node.inputs.pet_pvc_mni = []
            self.output_node.inputs.pet_pvc_suvr = []
            self.output_node.inputs.pet_pvc_suvr_masked = []
            self.output_node.inputs.pvc_atlas_statistics = []
