# coding: utf8

import clinica.pipelines.engine as cpe

__author__ = ["Alexandre Routier", "Thomas Jacquemont"]
__copyright__ = "Copyright 2016-2019 The Aramis Lab Team"
__credits__ = ["Nipype"]
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__status__ = "Development"

# Use hash instead of parameters for iterables folder names
# Otherwise path will be too long and generate OSError
from nipype import config
cfg = dict(execution={'parameterize_dirs': False})
config.update_config(cfg)


class DwiDti(cpe.Pipeline):
    """DTI-based processing of DWI datasets.

    Args:
        input_dir(str): Input directory in a CAPS hierarchy.
        output_dir(str): Output directory in a CAPS hierarchy.
        subjects_sessions_list(str): The Subjects-Sessions list file (in .tsv
            format).

    Returns:
        A clinica pipeline object containing the DwiDti pipeline.

    Raises:

    """

    def check_custom_dependencies(self): pass

    def get_input_fields(self):
        """Specify the list of possible inputs of this pipelines.

        Returns:
            A list of (string) input fields name.
        """
        input_list = ['preproc_dwi', 'preproc_bvec', 'preproc_bval', 'b0_mask']
        return input_list

    def get_output_fields(self):
        """Specify the list of possible outputs of this pipelines.

        Returns:
            A list of (string) output fields name.
        """
        output_list_dti = [
            'dti', 'fa', 'md', 'ad', 'rd', 'decfa',
            'registered_fa', 'registered_md', 'registered_ad', 'registered_rd',
            'statistics_fa', 'statistics_md', 'statistics_ad', 'statistics_rd',
            'b_spline_transform', 'affine_matrix'
        ]

        return output_list_dti

    def build_input_node(self):
        """Build and connect an input node to the pipelines.
        """
        import re
        from colorama import Fore
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        from clinica.iotools.grabcaps import CAPSLayout
        from clinica.utils.stream import cprint
        from clinica.utils.io import check_input_caps_files
        from clinica.utils.exceptions import ClinicaCAPSError

        # Remove 'sub-' prefix from participant IDs
        participants_regex = '|'.join(sub[4:] for sub in self.subjects)
        # Remove 'ses-' prefix from session IDs
        sessions_regex = '|'.join(ses[4:] for ses in self.sessions)

        caps_layout = CAPSLayout(self.caps_directory)

        error_message = ""
        # Inputs from dwi-preprocessing pipeline
        # ======================================
        # Preprocessed DWI
        dwi_files = caps_layout.get(type='dwi', suffix='preproc', extensions=['.nii|.nii.gz'], return_type='file',
                                    subject=participants_regex, session=sessions_regex)
        error_message += check_input_caps_files(dwi_files, "DWI_PREPROC_NII", "dwi-preprocessing",
                                                self.caps_directory, self.subjects, self.sessions)

        # B0 brainmask
        dwi_brainmask_files = caps_layout.get(type='dwi', suffix='brainmask', extensions=['.nii|.nii.gz'], return_type='file',
                                              subject=participants_regex, session=sessions_regex)
        error_message += check_input_caps_files(dwi_brainmask_files, "DWI_PREPROC_BM", "dwi-preprocessing",
                                                self.caps_directory, self.subjects, self.sessions)

        # Preprocessed bvec
        bvec_files = caps_layout.get(type='dwi', suffix='preproc', return_type='file', extensions='bvec',
                                     subject=participants_regex, session=sessions_regex)
        error_message += check_input_caps_files(bvec_files, "DWI_PREPROC_BVEC", "dwi-preprocessing",
                                                self.caps_directory, self.subjects, self.sessions)

        # Preprocessed bval
        bval_files = caps_layout.get(type='dwi', suffix='preproc', extensions='bval', return_type='file',
                                     subject=participants_regex, session=sessions_regex)
        error_message += check_input_caps_files(bval_files, "DWI_PREPROC_BVAL", "dwi-preprocessing",
                                                self.caps_directory, self.subjects, self.sessions)

        if error_message:
            raise ClinicaCAPSError(error_message)

        if len(dwi_files) == 0:
            import sys
            cprint('%s\nEither all the images were already run by the pipeline or no image was found to run the pipeline. '
                   'The program will now exit.%s' % (Fore.BLUE, Fore.RESET))
            sys.exit(0)
        else:
            p_id_images_to_process = [re.search(r'(sub-[a-zA-Z0-9]+)', caps_file).group() for caps_file in dwi_files]
            s_id_images_to_process = [re.search(r'(ses-[a-zA-Z0-9]+)', caps_file).group() for caps_file in dwi_files]
            images_to_process = ', '.join(p_id[4:] + '|' + s_id[4:]
                                          for p_id, s_id in zip(p_id_images_to_process, s_id_images_to_process))
            cprint('The pipeline will be run on the following subject(s): %s' % images_to_process)
            cprint('The pipeline will last approximately 20 minutes per image.')

        read_node = npe.Node(name="ReadingFiles",
                             iterables=[
                                 ('b0_mask', dwi_brainmask_files),
                                 ('preproc_dwi', dwi_files),
                                 ('preproc_bval', bval_files),
                                 ('preproc_bvec', bvec_files),
                             ],
                             synchronize=True,
                             interface=nutil.IdentityInterface(
                                     fields=self.get_input_fields()))
        self.connect([
            (read_node, self.input_node, [('b0_mask', 'b0_mask')]),
            (read_node, self.input_node, [('preproc_dwi', 'preproc_dwi')]),
            (read_node, self.input_node, [('preproc_bval', 'preproc_bval')]),
            (read_node, self.input_node, [('preproc_bvec', 'preproc_bvec')]),
        ])

    def build_output_node(self):
        """Build and connect an output node to the pipelines.
        """
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        import nipype.interfaces.io as nio
        from clinica.utils.io import fix_join

        import clinica.pipelines.dwi_dti.dwi_dti_utils as utils

        # Find container path from filename
        container_path = npe.Node(nutil.Function(
            input_names=['dwi_filename'],
            output_names=['container'],
            function=utils.dwi_container_from_filename),
            name='container_path')

        rename_into_caps = npe.Node(nutil.Function(
            input_names=['in_caps_dwi', 'in_norm_fa', 'in_norm_md',
                         'in_norm_ad', 'in_norm_rd', 'in_b_spline_transform',
                         'in_affine_matrix'],
            output_names=['out_caps_fa', 'out_caps_md', 'out_caps_ad',
                          'out_caps_rd', 'out_caps_b_spline_transform',
                          'out_caps_affine_matrix'],
            function=utils.rename_into_caps),
            name='rename_into_caps')

        # Writing results into CAPS
        write_results = npe.Node(name='write_results',
                                 interface=nio.DataSink())
        write_results.inputs.base_directory = self.caps_directory
        write_results.inputs.parameterization = False

        self.connect([
           (self.input_node, container_path, [('preproc_dwi', 'dwi_filename')]),

           (container_path,   write_results, [(('container', fix_join, 'dwi', 'dti_based_processing'), 'container')]),
           (self.output_node, write_results, [('dti',   'native_space.@dti')]),
           (self.output_node, write_results, [('fa',    'native_space.@fa'),
                                              ('md',    'native_space.@md'),
                                              ('ad',    'native_space.@ad'),
                                              ('rd',    'native_space.@rd'),
                                              ('decfa', 'native_space.@decfa')]),

           (self.input_node,  rename_into_caps, [('preproc_dwi',        'in_caps_dwi')]),
           (self.output_node, rename_into_caps, [('registered_fa',      'in_norm_fa'),
                                                 ('registered_md',      'in_norm_md'),
                                                 ('registered_ad',      'in_norm_ad'),
                                                 ('registered_rd',      'in_norm_rd'),
                                                 ('affine_matrix',      'in_affine_matrix'),
                                                 ('b_spline_transform', 'in_b_spline_transform')]),

           (rename_into_caps, write_results, [('out_caps_fa',                 'normalized_space.@registered_fa'),
                                              ('out_caps_md',                 'normalized_space.@registered_md'),
                                              ('out_caps_ad',                 'normalized_space.@registered_ad'),
                                              ('out_caps_rd',                 'normalized_space.@registered_rd'),
                                              ('out_caps_affine_matrix',      'normalized_space.@affine_matrix'),
                                              ('out_caps_b_spline_transform', 'normalized_space.@b_spline_transform')]),

           (self.output_node, write_results, [('statistics_fa', 'atlas_statistics.@statistics_fa'),
                                              ('statistics_md', 'atlas_statistics.@statistics_md'),
                                              ('statistics_ad', 'atlas_statistics.@statistics_ad'),
                                              ('statistics_rd', 'atlas_statistics.@statistics_rd')]) 
        ])

    def build_core_nodes(self):
        """Build and connect the core nodes of the pipelines.
        """
        import clinica.pipelines.dwi_dti.dwi_dti_workflows as workflows
        import clinica.pipelines.dwi_dti.dwi_dti_utils as utils

        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        import nipype.interfaces.fsl as fsl
        import nipype.interfaces.mrtrix as mrtrix
        from clinica.lib.nipype.interfaces.mrtrix3.utils import TensorMetrics
        from clinica.lib.nipype.interfaces.mrtrix.preprocess import DWI2Tensor

        # Nodes creation
        # ==============
        get_bids_identifier = npe.Node(interface=nutil.Function(
            input_names=['caps_dwi_filename'],
            output_names=['bids_identifier'],
            function=utils.extract_bids_identifier_from_caps_filename),
            name='0-Get_BIDS_Identifier')

        get_caps_filenames = npe.Node(interface=nutil.Function(
            input_names=['caps_dwi_filename'],
            output_names=['bids_source', 'out_dti',
                          'out_fa', 'out_md', 'out_ad', 'out_rd'],
            function=utils.get_caps_filenames),
            name='0-CAPS_Filenames')

        convert_gradients = npe.Node(
            interface=mrtrix.FSL2MRTrix(),
            name='0-Convert_FSL_Gradient')

        dwi_to_dti = npe.Node(
            # interface=mrtrix.DWI2Tensor(),
            interface=DWI2Tensor(),
            name='1-Compute_DTI')

        dti_to_metrics = npe.Node(
            # interface=mrtrix3.TensorMetrics(),
            interface=TensorMetrics(),
            name='2-DTI-based_Metrics')

        register_on_jhu_atlas = workflows.register_dti_maps_on_atlas(
            working_directory=self.base_dir,
            name="3-Register_DTI_Maps_On_JHU")

        scalar_analysis = npe.Node(
            interface=nutil.Function(
                input_names=['in_registered_map', 'name_map', 'prefix_file'],
                output_names=['atlas_statistics_list'],
                function=utils.statistics_on_atlases),
            name='4-Scalar_Analysis')
        scalar_analysis_fa = scalar_analysis.clone('4-Scalar_Analysis_FA')
        scalar_analysis_fa.inputs.name_map = 'FA'
        scalar_analysis_md = scalar_analysis.clone('4-Scalar_Analysis_MD')
        scalar_analysis_md.inputs.name_map = 'MD'
        scalar_analysis_ad = scalar_analysis.clone('4-Scalar_Analysis_AD')
        scalar_analysis_ad.inputs.name_map = 'AD'
        scalar_analysis_rd = scalar_analysis.clone('4-Scalar_Analysis_RD')
        scalar_analysis_rd.inputs.name_map = 'RD'

        thres_map = npe.Node(fsl.Threshold(thresh=0.0),
                             iterfield=['in_file'],
                             name='5-Remove_Negative')
        thres_fa = thres_map.clone('5-Remove_Negative_FA')
        thres_md = thres_map.clone('5-Remove_Negative_MD')
        thres_ad = thres_map.clone('5-Remove_Negative_AD')
        thres_rd = thres_map.clone('5-Remove_Negative_RD')

        print_begin_message = npe.Node(
            interface=nutil.Function(
                input_names=['in_bids_or_caps_file'],
                function=utils.print_begin_pipeline),
            name='Write-Begin_Message')

        print_end_message = npe.Node(
            interface=nutil.Function(
                input_names=['in_bids_or_caps_file', 'final_file_1', 'final_file_2'],
                function=utils.print_end_pipeline),
            name='Write-End_Message')

        # Connection
        # ==========
        self.connect([
            (self.input_node, get_caps_filenames, [('preproc_dwi', 'caps_dwi_filename')]),
            # Print begin message
            (self.input_node, print_begin_message, [('preproc_dwi', 'in_bids_or_caps_file')]),
            # Get BIDS/CAPS identifier from filename
            (self.input_node, get_bids_identifier, [('preproc_dwi', 'caps_dwi_filename')]),
            # Convert FSL gradient files (bval/bvec) to MRtrix format
            (self.input_node, convert_gradients, [('preproc_bval', 'bval_file'),
                                                  ('preproc_bvec', 'bvec_file')]),
            # Computation of the DTI model
            (self.input_node,    dwi_to_dti, [('b0_mask',       'in_mask'),
                                              ('preproc_dwi',   'in_file')]),
            (convert_gradients,  dwi_to_dti, [('encoding_file', 'encoding_file')]),
            (get_caps_filenames, dwi_to_dti, [('out_dti',       'out_filename')]),
            # Computation of the different metrics from the DTI
            (get_caps_filenames, dti_to_metrics, [('out_fa',  'out_fa')]),
            (get_caps_filenames, dti_to_metrics, [('out_md',  'out_adc')]),
            (get_caps_filenames, dti_to_metrics, [('out_ad',  'out_ad')]),
            (get_caps_filenames, dti_to_metrics, [('out_rd',  'out_rd')]),
            (self.input_node,    dti_to_metrics, [('b0_mask', 'in_mask')]),
            (dwi_to_dti,         dti_to_metrics, [('tensor',  'in_file')]),
            # Register DTI maps on JHU atlas
            (dti_to_metrics,    register_on_jhu_atlas, [('out_fa',  'inputnode.in_fa'),
                                                        ('out_adc', 'inputnode.in_md'),
                                                        ('out_ad',  'inputnode.in_ad'),
                                                        ('out_rd',  'inputnode.in_rd')]),
            # Generate regional TSV files
            (get_bids_identifier,   scalar_analysis_fa, [('bids_identifier',        'prefix_file')]),
            (register_on_jhu_atlas, scalar_analysis_fa, [('outputnode.out_norm_fa', 'in_registered_map')]),
            (get_bids_identifier,   scalar_analysis_md, [('bids_identifier',        'prefix_file')]),
            (register_on_jhu_atlas, scalar_analysis_md, [('outputnode.out_norm_md', 'in_registered_map')]),
            (get_bids_identifier,   scalar_analysis_ad, [('bids_identifier',        'prefix_file')]),
            (register_on_jhu_atlas, scalar_analysis_ad, [('outputnode.out_norm_ad', 'in_registered_map')]),
            (get_bids_identifier,   scalar_analysis_rd, [('bids_identifier',        'prefix_file')]),
            (register_on_jhu_atlas, scalar_analysis_rd, [('outputnode.out_norm_rd', 'in_registered_map')]),
            # Remove negative values from the DTI maps:
            (get_caps_filenames, thres_fa, [('out_fa',  'out_file')]),
            (dti_to_metrics,     thres_fa, [('out_fa',  'in_file')]),

            (get_caps_filenames, thres_md, [('out_md',  'out_file')]),
            (dti_to_metrics,     thres_md, [('out_adc', 'in_file')]),

            (get_caps_filenames, thres_ad, [('out_ad',  'out_file')]),
            (dti_to_metrics,     thres_ad, [('out_ad',  'in_file')]),

            (get_caps_filenames, thres_rd, [('out_rd',  'out_file')]),
            (dti_to_metrics,     thres_rd, [('out_rd',  'in_file')]),
            # Output node
            (dwi_to_dti,            self.output_node, [('tensor',   'dti')]),
            (thres_fa,              self.output_node, [('out_file', 'fa')]),
            (thres_md,              self.output_node, [('out_file', 'md')]),
            (thres_ad,              self.output_node, [('out_file', 'ad')]),
            (thres_rd,              self.output_node, [('out_file', 'rd')]),
            (dti_to_metrics,        self.output_node, [('out_evec', 'decfa')]),
            (register_on_jhu_atlas, self.output_node, [('outputnode.out_norm_fa',            'registered_fa'),
                                                       ('outputnode.out_norm_md',            'registered_md'),
                                                       ('outputnode.out_norm_ad',            'registered_ad'),
                                                       ('outputnode.out_norm_rd',            'registered_rd'),
                                                       ('outputnode.out_affine_matrix',      'affine_matrix'),
                                                       ('outputnode.out_b_spline_transform', 'b_spline_transform')]),
            (scalar_analysis_fa,    self.output_node, [('atlas_statistics_list', 'statistics_fa')]),
            (scalar_analysis_md,    self.output_node, [('atlas_statistics_list', 'statistics_md')]),
            (scalar_analysis_ad,    self.output_node, [('atlas_statistics_list', 'statistics_ad')]),
            (scalar_analysis_rd,    self.output_node, [('atlas_statistics_list', 'statistics_rd')]), 
            # Print end message
            (self.input_node,    print_end_message, [('preproc_dwi',           'in_bids_or_caps_file')]),
            (thres_rd,           print_end_message, [('out_file',              'final_file_1')]),
            (scalar_analysis_rd, print_end_message, [('atlas_statistics_list', 'final_file_2')]),
        ])
