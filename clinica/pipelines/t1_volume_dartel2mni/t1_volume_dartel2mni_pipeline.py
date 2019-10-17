# coding: utf8

import clinica.pipelines.engine as cpe

__author__ = "Jorge Samper-Gonzalez"
__copyright__ = "Copyright 2016-2019 The Aramis Lab Team"
__credits__ = ["Jorge Samper-Gonzalez"]
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__maintainer__ = "Jorge Samper-Gonzalez"
__email__ = "jorge.samper-gonzalez@inria.fr"
__status__ = "Development"


class T1VolumeDartel2MNI(cpe.Pipeline):
    """T1VolumeDartel2MNI - Dartel template to MNI.

    Args:
        input_dir: A BIDS directory.
        output_dir: An empty output directory where CAPS structured data will be written.
        subjects_sessions_list: The Subjects-Sessions list file (in .tsv format).

    Returns:
        A clinica pipeline object containing the T1VolumeDartel2MNI pipeline.
    """
    def __init__(self, bids_directory=None, caps_directory=None, tsv_file=None, name=None, group_id='default'):
        import os
        from clinica.utils.exceptions import ClinicaException

        super(T1VolumeDartel2MNI, self).__init__(bids_directory, caps_directory, tsv_file, name)

        if not group_id.isalnum():
            raise ValueError('Not valid group_id value. It must be composed only by letters and/or numbers')

        self._group_id = group_id
        # Check that group already exists
        if not os.path.exists(os.path.join(os.path.abspath(caps_directory), 'groups', 'group-' + group_id)):
            error_message = group_id \
                            + ' does not exists, please choose another one (or maybe you need to run t1-volume-create-dartel).' \
                            + '\nGroups that already exist in your CAPS directory are: \n'
            list_groups = os.listdir(os.path.join(os.path.abspath(caps_directory), 'groups'))
            is_empty = True
            for e in list_groups:
                if e.startswith('group-'):
                    error_message += e + ' \n'
                    is_empty = False
            if is_empty is True:
                error_message += 'NO GROUP FOUND'
            raise ClinicaException(error_message)

        # Default parameters
        self._parameters = {'tissues': [1, 2, 3],
                            'bounding_box': None,
                            'voxel_size': None,
                            'modulation': True,
                            'fwhm': [8]
                            # 'atlas_list': ['AAL2', 'LPBA40', 'Neuromorphometrics', 'AICHA', 'Hammers']
                            }

    def check_custom_dependencies(self):
        """Check dependencies that can not be listed in the `info.json` file.
        """
        pass

    def get_input_fields(self):
        """Specify the list of possible inputs of this pipelines.

        Returns:
            A list of (string) input fields name.
        """

        return ['native_segmentations', 'flowfield_files', 'template_file']

    def get_output_fields(self):
        """Specify the list of possible outputs of this pipelines.

        Returns:
            A list of (string) output fields name.
        """

        return ['normalized_files', 'smoothed_normalized_files', 'atlas_statistics']

    def build_input_node(self):
        """Build and connect an input node to the pipelines.
        """

        import nipype.pipeline.engine as npe
        import nipype.interfaces.io as nio
        import nipype.interfaces.utility as nutil
        from clinica.utils.inputs import clinica_file_reader, clinica_group_reader
        from clinica.utils.exceptions import ClinicaCAPSError, ClinicaException
        from colorama import Fore

        tissue_names = {1: 'graymatter',
                        2: 'whitematter',
                        3: 'csf',
                        4: 'bone',
                        5: 'softtissue',
                        6: 'background'
                        }

        all_errors = []
        read_input_node = npe.Node(name="LoadingCLIArguments",
                                   interface=nutil.IdentityInterface(
                                       fields=self.get_input_fields(),
                                       mandatory_inputs=True))

        # Segmented Tissues
        # =================
        tissues_input = []
        for tissue_number in self.parameters['tissues']:
            try:
                current_file = clinica_file_reader(self.subjects,
                                                   self.sessions,
                                                   self.caps_directory,
                                                   {'pattern': 't1/spm/segmentation/native_space/*_*_T1w_segm-'
                                                               + tissue_names[tissue_number] + '_probability.nii*',
                                                    'description': 'SPM based probability of ' + tissue_names[tissue_number]
                                                                   + ' based on T1w-MRI in native space',
                                                    'needed_pipeline': 't1-volume-tissue-segmentation'})
                tissues_input.append(current_file)
            except ClinicaException as e:
                all_errors.append(e)
        # Tissues_input has a length of len(self.parameters['mask_tissues']). Each of these elements has a size of
        # len(self.subjects). We want the opposite : a list of size len(self.subjects) whose elements have a size of
        # len(self.parameters['mask_tissues']. The trick is to iter on elements with zip(*mylist)
        tissues_input_rearranged = []
        for subject_tissue_list in zip(*tissues_input):
            tissues_input_rearranged.append(subject_tissue_list)

        read_input_node.inputs.native_segmentations = tissues_input_rearranged

        # Flow Fields
        # ===========
        try:
            read_input_node.inputs.flowfield_files = clinica_file_reader(self.subjects,
                                                                         self.sessions,
                                                                         self.caps_directory,
                                                                         {'pattern': 't1/spm/dartel/group-' + self._group_id
                                                                                     + '/sub-*_ses-*_T1w_target-' + self._group_id
                                                                                     + '_transformation-forward_deformation.nii*',
                                                                          'description': 'flowfield files (forward transformation) from native space to '
                                                                                         + self._group_id + ' space',
                                                                          'needed_pipeline': 't1-volume-create-dartel'})
        except ClinicaException as e:
            all_errors.append(e)

        # Dartel Template
        # ================
        try:
            read_input_node.inputs.template_file = clinica_group_reader(self.caps_directory,
                                                                        {'pattern': 'group-' + self._group_id
                                                                                    + '/t1/group-' + self._group_id
                                                                                    + '_template.nii*',
                                                                         'description': 'T1w template file of group '
                                                                                        + self._group_id,
                                                                         'needed_pipeline': 't1-volume or t1-volume-create-dartel'})
        except ClinicaException as e:
            all_errors.append(e)

        if len(all_errors) > 0:
            error_message = 'Clinica faced error(s) while trying to read files in your CAPS/BIDS directories.\n'
            for msg in all_errors:
                error_message += str(msg)
            raise ClinicaCAPSError(error_message)

        self.connect([
            (read_input_node, self.input_node, [('native_segmentations', 'native_segmentations')]),
            (read_input_node, self.input_node, [('flowfield_files', 'flowfield_files')]),
            (read_input_node, self.input_node, [('template_file', 'template_file')])
        ])

    def build_output_node(self):
        """Build and connect an output node to the pipelines.
        """
        import os.path as op
        import nipype.pipeline.engine as npe
        import nipype.interfaces.io as nio
        import re
        from clinica.utils.io import zip_nii

        # Writing normalized images (and smoothed) into CAPS
        # ==================================================
        write_normalized_node = npe.MapNode(name='write_normalized_node',
                                            iterfield=['container', 'normalized_files', 'smoothed_normalized_files'],
                                            interface=nio.DataSink(infields=['normalized_files',
                                                                             'smoothed_normalized_files']))
        write_normalized_node.inputs.base_directory = self.caps_directory
        write_normalized_node.inputs.parameterization = False
        write_normalized_node.inputs.container = ['subjects/' + self.subjects[i] + '/' + self.sessions[i] +
                                                  '/t1/spm/dartel/group-' + self._group_id
                                                  for i in range(len(self.subjects))]
        write_normalized_node.inputs.regexp_substitutions = [
            (r'(.*)c1(sub-.*)(\.nii(\.gz)?)$', r'\1\2_segm-graymatter_probability\3'),
            (r'(.*)c2(sub-.*)(\.nii(\.gz)?)$', r'\1\2_segm-whitematter_probability\3'),
            (r'(.*)c3(sub-.*)(\.nii(\.gz)?)$', r'\1\2_segm-csf_probability\3'),
            (r'(.*)c4(sub-.*)(\.nii(\.gz)?)$', r'\1\2_segm-bone_probability\3'),
            (r'(.*)c5(sub-.*)(\.nii(\.gz)?)$', r'\1\2_segm-softtissue_probability\3'),
            (r'(.*)c6(sub-.*)(\.nii(\.gz)?)$', r'\1\2_segm-background_probability\3'),
            (r'(.*)mw(sub-.*)_probability(\.nii(\.gz)?)$', r'\1\2_space-Ixi549Space_modulated-on_probability\3'),
            (r'(.*)w(sub-.*)_probability(\.nii(\.gz)?)$', r'\1\2_space-Ixi549Space_modulated-off_probability\3'),
            (r'(.*)/normalized_files/(sub-.*)$', r'\1/\2'),
            (r'(.*)/smoothed_normalized_files/(fwhm-[0-9]+mm)_(sub-.*)_probability(\.nii(\.gz)?)$',
             r'\1/\3_\2_probability\4'),
            (r'trait_added', r'')
        ]
        self.connect([
            (self.output_node, write_normalized_node, [(('normalized_files', zip_nii, True), 'normalized_files'),
                                                       (('smoothed_normalized_files', zip_nii, True),
                                                        'smoothed_normalized_files')])
        ])

    def build_core_nodes(self):
        """Build and connect the core nodes of the pipelines.
        """
        import nipype.interfaces.spm as spm
        import nipype.pipeline.engine as npe
        import nipype.interfaces.utility as nutil
        from clinica.utils.io import unzip_nii
        from ..t1_volume_dartel2mni import t1_volume_dartel2mni_utils as dartel2mni_utils
        from clinica.utils.spm import get_tpm

        # Get Tissue Probability Map from SPM
        tissue_map = get_tpm()

        # Unzipping
        # =========
        unzip_tissues_node = npe.MapNode(nutil.Function(input_names=['in_file'],
                                                        output_names=['out_file'],
                                                        function=unzip_nii),
                                         name='unzip_tissues_node',
                                         iterfield=['in_file'])
        unzip_flowfields_node = npe.MapNode(nutil.Function(input_names=['in_file'],
                                                           output_names=['out_file'],
                                                           function=unzip_nii),
                                            name='unzip_flowfields_node',
                                            iterfield=['in_file'])
        unzip_template_node = npe.Node(nutil.Function(input_names=['in_file'],
                                                      output_names=['out_file'],
                                                      function=unzip_nii),
                                       name='unzip_template_node')

        # DARTEL2MNI Registration
        # =======================
        dartel2mni_node = npe.MapNode(spm.DARTELNorm2MNI(),
                                      name='dartel2MNI',
                                      iterfield=['apply_to_files', 'flowfield_files'])

        if self.parameters['bounding_box'] is not None:
            dartel2mni_node.inputs.bounding_box = self.parameters['bounding_box']
        if self.parameters['voxel_size'] is not None:
            dartel2mni_node.inputs.voxel_size = self.parameters['voxel_size']
        dartel2mni_node.inputs.modulate = self.parameters['modulation']
        dartel2mni_node.inputs.fwhm = 0

        # Smoothing
        # =========
        if self.parameters['fwhm'] is not None and len(self.parameters['fwhm']) > 0:
            smoothing_node = npe.MapNode(spm.Smooth(),
                                         name='smoothing_node',
                                         iterfield=['in_files'])

            smoothing_node.iterables = [('fwhm', [[x, x, x] for x in self.parameters['fwhm']]),
                                        ('out_prefix', ['fwhm-' + str(x) + 'mm_' for x in self.parameters['fwhm']])]
            smoothing_node.synchronize = True

            join_smoothing_node = npe.JoinNode(interface=nutil.Function(input_names=['smoothed_normalized_files'],
                                                                        output_names=['smoothed_normalized_files'],
                                                                        function=dartel2mni_utils.join_smoothed_files),
                                               joinsource='smoothing_node',
                                               joinfield='smoothed_normalized_files',
                                               name='join_smoothing_node')
            self.connect([
                (dartel2mni_node, smoothing_node, [('normalized_files', 'in_files')]),
                (smoothing_node, join_smoothing_node, [('smoothed_files', 'smoothed_normalized_files')]),
                (join_smoothing_node, self.output_node, [('smoothed_normalized_files', 'smoothed_normalized_files')])
            ])
        else:
            self.output_node.inputs.smoothed_normalized_files = []

        # Connection
        # ==========
        self.connect([
            (self.input_node, unzip_tissues_node, [('native_segmentations', 'in_file')]),
            (self.input_node, unzip_flowfields_node, [('flowfield_files', 'in_file')]),
            (self.input_node, unzip_template_node, [('template_file', 'in_file')]),
            (unzip_tissues_node, dartel2mni_node, [('out_file', 'apply_to_files')]),
            (unzip_flowfields_node, dartel2mni_node, [(('out_file', dartel2mni_utils.prepare_flowfields, self.parameters['tissues']),
                                                       'flowfield_files')]),
            (unzip_template_node, dartel2mni_node, [('out_file', 'template_file')]),
            (dartel2mni_node, self.output_node, [('normalized_files', 'normalized_files')])
        ])
