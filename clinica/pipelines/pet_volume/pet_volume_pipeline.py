# coding: utf8

import clinica.pipelines.engine as cpe
from nipype import config

# Use hash instead of parameters for iterables folder names
# Otherwise path will be too long and generate OSError
cfg = dict(execution={'parameterize_dirs': False})
config.update_config(cfg)


class PETVolume(cpe.Pipeline):
    """PETVolume - Volume-based processing of PET images using SPM.

    Returns:
        A clinica pipeline object containing the PETVolume pipeline.
    """
    def check_pipeline_parameters(self):
        """Check pipeline parameters."""
        from clinica.utils.group import check_group_label
        default_atlases = ['AAL2', 'LPBA40', 'Neuromorphometrics', 'AICHA', 'Hammers']

        if 'group_id' not in self.parameters.keys():
            raise KeyError('Missing compulsory group_id key in pipeline parameter.')
        if 'psf_tsv' not in self.parameters.keys():
            self.parameters['psf_tsv'] = None
        if 'pet_tracer' not in self.parameters.keys():
            self.parameters['pet_tracer'] = 'fdg'
        if 'mask_tissues' not in self.parameters.keys():
            self.parameters['mask_tissues'] = [1, 2, 3]
        if 'mask_threshold' not in self.parameters.keys():
            self.parameters['mask_threshold'] = 0.3
        if 'pvc_mask_tissues' not in self.parameters.keys():
            self.parameters['pvc_mask_tissues'] = [1, 2, 3]
        if 'smooth' not in self.parameters.keys():
            self.parameters['smooth'] = [8]
        if 'atlases' not in self.parameters.keys():
            self.parameters['atlases'] = default_atlases

        check_group_label(self.parameters['group_id'])

    def check_custom_dependencies(self):
        """Check dependencies that can not be listed in the `info.json` file.
        """
        pass

    def get_input_fields(self):
        """Specify the list of possible inputs of this pipelines.

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
        """Specify the list of possible outputs of this pipelines.

        Returns:
            A list of (string) output fields name.
        """

        return ['pet_t1_native',
                'pet_mni',
                'pet_suvr',
                'binary_mask',
                'pet_suvr_masked',
                'pet_suvr_masked_smoothed',
                'pet_pvc',
                'pet_pvc_mni',
                'pet_pvc_suvr',
                'pet_pvc_suvr_masked',
                'pet_pvc_suvr_masked_smoothed',
                'atlas_statistics',
                'pvc_atlas_statistics'
                ]

    def build_input_node(self):
        """Build and connect an input node to the pipelines.
        """
        import os
        from os.path import join, split, realpath, exists
        from colorama import Fore
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        from .pet_volume_utils import read_psf_information
        from clinica.utils.inputs import clinica_file_reader, clinica_group_reader, insensitive_glob
        from clinica.utils.input_files import (t1_volume_final_group_template,
                                               t1_volume_native_tpm,
                                               t1_volume_native_tpm_in_mni,
                                               t1_volume_deformation_to_template)
        from clinica.utils.exceptions import ClinicaException
        from clinica.utils.ux import print_groups_in_caps_directory, print_images_to_process
        from clinica.iotools.utils.data_handling import check_relative_volume_location_in_world_coordinate_system
        import clinica.utils.input_files as input_files
        from clinica.utils.filemanip import save_participants_sessions
        from clinica.utils.stream import cprint

        # Check that group already exists
        if not exists(join(self.caps_directory, 'groups', 'group-' + self.parameters['group_id'])):
            print_groups_in_caps_directory(self.caps_directory)
            raise ClinicaException(
                '%sGroup %s does not exist. Did you run t1-volume or t1-volume-create-dartel pipeline?%s' %
                (Fore.RED, self.parameters['group_id'], Fore.RESET)
            )

        # Tissues DataGrabber
        # ====================
        all_errors = []

        # Grab reference mask
        if self.parameters['pet_tracer'] == 'fdg':
            reference_mask_template = 'region-pons_eroded-6mm_mask.nii*'
            self.parameters['suvr_region'] = 'pons'
            pet_file_to_grab = input_files.PET_FDG_NII
        elif self.parameters['pet_tracer'] == 'av45':
            reference_mask_template = 'region-cerebellumPons_eroded-6mm_mask.nii*'
            self.parameters['suvr_region'] = 'cerebellumPons'
            pet_file_to_grab = input_files.PET_AV45_NII
        else:
            raise NotImplementedError('Only "fdg" or "av45" tracers are currently accepted (given tracer "%s").' %
                                      self.parameters['pet_tracer'])
        reference_mask_path = join(split(realpath(__file__))[0], '../../resources/masks', reference_mask_template)
        reference_mask_file = insensitive_glob(reference_mask_path)
        if len(reference_mask_file) != 1:
            all_errors.append('An error happened while trying to get reference mask ' + str(reference_mask_path))
        else:
            reference_mask_file = reference_mask_file[0]

        # PET from BIDS directory
        try:
            pet_bids = clinica_file_reader(self.subjects,
                                           self.sessions,
                                           self.bids_directory,
                                           pet_file_to_grab)
        except ClinicaException as e:
            all_errors.append(e)

        # Native T1w-MRI
        try:
            t1w_bids = clinica_file_reader(self.subjects,
                                           self.sessions,
                                           self.bids_directory,
                                           input_files.T1W_NII)
        except ClinicaException as e:
            all_errors.append(e)

        # mask_tissues
        tissues_input = []
        for tissue_number in self.parameters['mask_tissues']:
            try:
                current_file = clinica_file_reader(self.subjects,
                                                   self.sessions,
                                                   self.caps_directory,
                                                   t1_volume_native_tpm_in_mni(tissue_number, False))
                tissues_input.append(current_file)
            except ClinicaException as e:
                all_errors.append(e)
        # Tissues_input has a length of len(self.parameters['mask_tissues']). Each of these elements has a size of
        # len(self.subjects). We want the opposite: a list of size len(self.subjects) whose elements have a size of
        # len(self.parameters['mask_tissues']. The trick is to iter on elements with zip(*my_list)
        tissues_input_final = []
        for subject_tissue_list in zip(*tissues_input):
            tissues_input_final.append(subject_tissue_list)
        tissues_input = tissues_input_final

        # Flowfields
        try:
            flowfields_caps = clinica_file_reader(self.subjects,
                                                  self.sessions,
                                                  self.caps_directory,
                                                  t1_volume_deformation_to_template(self.parameters['group_id']))
        except ClinicaException as e:
            all_errors.append(e)

        # Dartel Template
        try:
            final_template = clinica_group_reader(self.caps_directory,
                                                  t1_volume_final_group_template(self.parameters['group_id']))
        except ClinicaException as e:
            all_errors.append(e)

        if self.parameters['psf_tsv'] is not None:
            iterables_fwhm = read_psf_information(self.parameters['psf_tsv'], self.subjects, self.sessions)
            self.parameters['apply_pvc'] = True
        else:
            iterables_fwhm = [[]] * len(self.subjects)
            self.parameters['apply_pvc'] = False

        if self.parameters['apply_pvc']:
            # pvc tissues input
            pvc_tissues_input = []
            for tissue_number in self.parameters['pvc_mask_tissues']:
                try:
                    current_file = clinica_file_reader(self.subjects,
                                                       self.sessions,
                                                       self.caps_directory,
                                                       t1_volume_native_tpm(tissue_number))
                    pvc_tissues_input.append(current_file)
                except ClinicaException as e:
                    all_errors.append(e)

            if len(all_errors) == 0:
                pvc_tissues_input_final = []
                for subject_tissue_list in zip(*pvc_tissues_input):
                    pvc_tissues_input_final.append(subject_tissue_list)
                pvc_tissues_input = pvc_tissues_input_final
        else:
            pvc_tissues_input = []

        if len(all_errors) > 0:
            error_message = 'Clinica faced error(s) while trying to read files in your CAPS/BIDS directories.\n'
            for msg in all_errors:
                error_message += str(msg)
            raise ClinicaException(error_message)

        check_relative_volume_location_in_world_coordinate_system('T1w-MRI', t1w_bids,
                                                                  self.parameters['pet_tracer'] + ' PET', pet_bids,
                                                                  self.bids_directory,
                                                                  self.parameters['pet_tracer'])

        # Save subjects to process in <WD>/<Pipeline.name>/participants.tsv
        folder_participants_tsv = os.path.join(self.base_dir, self.name)
        save_participants_sessions(self.subjects, self.sessions, folder_participants_tsv)

        if len(self.subjects):
            print_images_to_process(self.subjects, self.sessions)
            cprint('List available in %s' % os.path.join(folder_participants_tsv, 'participants.tsv'))
            cprint('The pipeline will last approximately 10 minutes per image.')

        read_input_node = npe.Node(name="LoadingCLIArguments",
                                   interface=nutil.IdentityInterface(
                                       fields=self.get_input_fields(),
                                       mandatory_inputs=True),
                                   iterables=[('pet_image', pet_bids),
                                              ('t1_image_native', t1w_bids),
                                              ('mask_tissues', tissues_input),
                                              ('fwhm', iterables_fwhm),
                                              ('flow_fields', flowfields_caps),
                                              ('pvc_mask_tissues', pvc_tissues_input)],
                                   synchronize=True)

        read_input_node.inputs.reference_mask = reference_mask_file
        read_input_node.inputs.dartel_template = final_template

        self.connect([(read_input_node, self.input_node, [('pet_image', 'pet_image')]),
                      (read_input_node, self.input_node, [('t1_image_native', 't1_image_native')]),
                      (read_input_node, self.input_node, [('mask_tissues', 'mask_tissues')]),
                      (read_input_node, self.input_node, [('flow_fields', 'flow_fields')]),
                      (read_input_node, self.input_node, [('dartel_template', 'dartel_template')]),
                      (read_input_node, self.input_node, [('reference_mask', 'reference_mask')]),
                      (read_input_node, self.input_node, [('fwhm', 'fwhm')]),
                      (read_input_node, self.input_node, [('pvc_mask_tissues', 'pvc_mask_tissues')])])

    def build_output_node(self):
        """Build and connect an output node to the pipelines.
        """

        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        import nipype.interfaces.io as nio
        from clinica.utils.filemanip import zip_nii
        from clinica.utils.nipype import fix_join
        import re
        import clinica.pipelines.pet_volume.pet_volume_utils as utils

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
            (r'(.*/)pet_suvr/suvr_wr(sub-.*)(\.nii(\.gz)?)$',
             r'\1\2_space-Ixi549Space_suvr-' + re.escape(self.parameters['suvr_region']) + r'_pet\3'),
            (r'(.*/)pet_pvc_suvr/suvr_wpvc-rbv_r(sub-.*)(\.nii(\.gz)?)$',
             r'\1\2_space-Ixi549Space_pvc-rbv_suvr-' + re.escape(self.parameters['suvr_region']) + r'_pet\3'),
            (r'(.*/)pet_suvr_masked/masked_suvr_wr(sub-.*)(\.nii(\.gz)?)$',
             r'\1\2_space-Ixi549Space_suvr-' + re.escape(self.parameters['suvr_region']) + r'_mask-brain_pet\3'),
            (r'(.*/)pet_pvc_suvr_masked/masked_suvr_wpvc-rbv_r(sub-.*)(\.nii(\.gz)?)$',
             r'\1\2_space-Ixi549Space_pvc-rbv_suvr-' + re.escape(self.parameters['suvr_region']) + r'_mask-brain_pet\3'),
            (r'(.*/)pet_suvr_masked_smoothed/(fwhm-[0-9]+mm)_masked_suvr_wr(sub-.*)(\.nii(\.gz)?)$',
             r'\1\3_space-Ixi549Space_suvr-' + re.escape(self.parameters['suvr_region']) + r'_mask-brain_\2_pet\4'),
            (r'(.*/)pet_pvc_suvr_masked_smoothed/(fwhm-[0-9]+mm)_masked_suvr_wpvc-rbv_r(sub-.*)(\.nii(\.gz)?)$',
             r'\1\3_space-Ixi549Space_pvc-rbv_suvr-' + re.escape(self.parameters['suvr_region']) + r'_mask-brain_\2_pet\4'),
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
             re.escape(self.parameters['suvr_region']) + r'\3'),
            (r'(.*/)pvc_(atlas_statistics/)suvr_wpvc-rbv_r(sub-.*)(_statistics\.tsv)$', r'\1\2\3' + r'_pvc-rbv_suvr-' +
             re.escape(self.parameters['suvr_region']) + r'\4')
        ]

        self.connect([(self.input_node, container_path, [('pet_image', 'pet_filename')]),
                      (container_path, write_images_node, [(('container', fix_join, 'group-' + self.parameters['group_id']),
                                                            'container')]),
                      (self.output_node, write_images_node, [(('pet_t1_native', zip_nii, True), 'pet_t1_native'),
                                                             (('pet_mni', zip_nii, True), 'pet_mni'),
                                                             (('pet_suvr', zip_nii, True), 'pet_suvr'),
                                                             (('binary_mask', zip_nii, True), 'binary_mask'),
                                                             (('pet_suvr_masked', zip_nii, True), 'pet_suvr_masked'),
                                                             (('pet_suvr_masked_smoothed', zip_nii, True),
                                                              'pet_suvr_masked_smoothed'),
                                                             (('pet_pvc', zip_nii, True), 'pet_pvc'),
                                                             (('pet_pvc_mni', zip_nii, True), 'pet_pvc_mni'),
                                                             (('pet_pvc_suvr', zip_nii, True), 'pet_pvc_suvr'),
                                                             (('pet_pvc_suvr_masked', zip_nii, True),
                                                              'pet_pvc_suvr_masked'),
                                                             (('pet_pvc_suvr_masked_smoothed', zip_nii, True),
                                                              'pet_pvc_suvr_masked_smoothed')]),
                      (container_path, write_atlas_node, [(('container', fix_join, 'group-' + self.parameters['group_id']),
                                                           'container')]),
                      (self.output_node, write_atlas_node, [('atlas_statistics', 'atlas_statistics'),
                                                            ('pvc_atlas_statistics', 'pvc_atlas_statistics')])
                      ])

    def build_core_nodes(self):
        """Build and connect the core nodes of the pipelines.
        """

        import nipype.interfaces.spm as spm
        import nipype.interfaces.spm.utils as spmutils
        from nipype.interfaces.petpvc import PETPVC
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        from clinica.utils.filemanip import unzip_nii
        import clinica.pipelines.pet_volume.pet_volume_utils as utils

        import os

        if 'SPMSTANDALONE_HOME' in os.environ:
            if 'MCR_HOME' in os.environ:
                matlab_cmd = (
                        os.path.join(
                            os.environ['SPMSTANDALONE_HOME'], 'run_spm12.sh')
                        + ' ' + os.environ['MCR_HOME']
                        + ' script')
                spm.SPMCommand.set_mlab_paths(matlab_cmd=matlab_cmd, use_mcr=True)

        # Initialize pipeline
        # ===================
        init_node = npe.Node(interface=nutil.Function(
            input_names=['pet_nii'],
            output_names=['pet_nii'],
            function=utils.init_input_node),
            name='init_pipeline')

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

        # Smoothing
        # =========
        if self.parameters['smooth'] is not None and len(self.parameters['smooth']) > 0:
            smoothing_node = npe.MapNode(spm.Smooth(),
                                         name='smoothing_node',
                                         iterfield=['fwhm', 'out_prefix'])
            smoothing_node.inputs.fwhm = [[x, x, x] for x in self.parameters['smooth']]
            smoothing_node.inputs.out_prefix = ['fwhm-' + str(x) + 'mm_' for x in self.parameters['smooth']]
            self.connect([
                (apply_mask, smoothing_node, [('masked_image_path', 'in_files')]),
                (smoothing_node, self.output_node, [('smoothed_files', 'pet_suvr_masked_smoothed')])
            ])
        else:
            self.output_node.inputs.pet_suvr_masked_smoothed = [[]]

        # Atlas Statistics
        # ================
        atlas_stats_node = npe.MapNode(nutil.Function(input_names=['in_image',
                                                                   'in_atlas_list'],
                                                      output_names=['atlas_statistics'],
                                                      function=utils.atlas_statistics),
                                       name='atlas_stats_node',
                                       iterfield=['in_image'])
        atlas_stats_node.inputs.in_atlas_list = self.parameters['atlases']

        # Connection
        # ==========
        self.connect([(self.input_node, init_node, [('pet_image', 'pet_nii')]),
                      (init_node, unzip_pet_image, [('pet_nii', 'in_file')]),
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
        if self.parameters['apply_pvc']:
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
            # Smoothing
            # =========
            if self.parameters['smooth'] is not None and len(self.parameters['smooth']) > 0:
                smoothing_pvc = npe.MapNode(spm.Smooth(),
                                            name='smoothing_pvc',
                                            iterfield=['fwhm', 'out_prefix'])
                smoothing_pvc.inputs.fwhm = [[x, x, x] for x in self.parameters['smooth']]
                smoothing_pvc.inputs.out_prefix = ['fwhm-' + str(x) + 'mm_' for x in self.parameters['smooth']]
                self.connect([
                    (apply_mask_pvc, smoothing_pvc, [('masked_image_path', 'in_files')]),
                    (smoothing_pvc, self.output_node, [('smoothed_files', 'pet_pvc_suvr_masked_smoothed')])
                ])
            else:
                self.output_node.inputs.pet_pvc_suvr_masked_smoothed = [[]]
            # Atlas Statistics
            # ================
            atlas_stats_pvc = npe.MapNode(nutil.Function(input_names=['in_image',
                                                                      'in_atlas_list'],
                                                         output_names=['atlas_statistics'],
                                                         function=utils.atlas_statistics),
                                          name='atlas_stats_pvc',
                                          iterfield=['in_image'])
            atlas_stats_pvc.inputs.in_atlas_list = self.parameters['atlases']

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
            self.output_node.inputs.pet_pvc = [[]]
            self.output_node.inputs.pet_pvc_mni = [[]]
            self.output_node.inputs.pet_pvc_suvr = [[]]
            self.output_node.inputs.pet_pvc_suvr_masked = [[]]
            self.output_node.inputs.pvc_atlas_statistics = [[]]
            self.output_node.inputs.pet_pvc_suvr_masked_smoothed = [[]]
