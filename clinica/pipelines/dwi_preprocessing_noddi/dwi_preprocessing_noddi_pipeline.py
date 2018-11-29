# coding: utf8

__author__ = "Junhao Wen"
__copyright__ = "Copyright 2016-2018, The Aramis Lab Team"
__credits__ = ["Junhao Wen"]
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__maintainer__ = "Junhao Wen"
__email__ = "Junhao.Wen@inria.fr"
__status__ = "Development"

# command line tool.
import clinica.pipelines.engine as cpe


class DwiPreprocessingNoddi(cpe.Pipeline):
    """dwi_preprocessing_noddi SHORT DESCRIPTION.

    Warnings:
        - A WARNING.

    Todos:
        - [x] A FILLED TODO ITEM.
        - [ ] AN ON-GOING TODO ITEM.

    Args: input_dir: A BIDS directory.  output_dir: An empty output directory
    where CAPS structured data will be written.  subjects_sessions_list: The
    Subjects-Sessions list file (in .tsv format).

    Returns:
        A clinica pipeline object containing the dwi_preprocessing_noddi
        pipeline.

    Raises:


    Example:
        >>> pipeline = DwiPreprocessingNoddi('~/MYDATASET_BIDS', '~/MYDATASET_CAPS')
        >>> pipeline.parameters = {
        >>>     # ...
        >>> }
        >>> pipeline.base_dir = '/tmp/'
        >>> pipeline.run()
    """


    def check_custom_dependencies(self):
        """Check dependencies that can not be listed in the `info.json` file.
        """
        pass


    def get_input_fields(self):
        """Specify the list of possible inputs of this pipeline.

        Returns:
            A list of (string) input fields name.
        """

        return ['bids_ap_dwi', 
                'bids_ap_dwi_bvec', 
                'bids_ap_dwi_bval', 
                'bids_pa_dwi', 
                'bids_pa_dwi_bvec', 
                'bids_pa_dwi_bval',
                'epi_param', 
                'epi_param_alt']  # Fill here the list


    def get_output_fields(self):
        """Specify the list of possible outputs of this pipeline.

        Returns:
            A list of (string) output fields name.
        """

        return ['preproc_dwi', 'preproc_bvec', 'preproc_bval', 'b0_mask']  # Fill here the list


    def build_input_node(self):
        """Build and connect an input node to the pipeline.
        """

        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        import dwi_preprocessing_noddi_utils as utils

        bids_ap_dwi, bids_ap_dwi_bvec, bids_ap_dwi_bval, bids_pa_dwi, bids_pa_dwi_bvec, bids_pa_dwi_bval = utils.grab_noddi_bids_files(self.bids_directory, self.tsv_file)
        list_tuple = utils.create_list_tuple(bids_ap_dwi, bids_ap_dwi_bvec, bids_ap_dwi_bval, bids_pa_dwi, bids_pa_dwi_bvec, bids_pa_dwi_bval)
        parallelsubjects = npe.Node(name="noddi_preprocessing_parallelization", interface=nutil.IdentityInterface(
            fields=['bids_ap_dwi', 'bids_ap_dwi_bvec', 'bids_ap_dwi_bval', 'bids_pa_dwi', 'bids_pa_dwi_bvec',
                    'bids_pa_dwi_bval'], mandatory_inputs=True), synchronize=True, iterables=[('bids_ap_dwi', 'bids_ap_dwi_bvec', 'bids_ap_dwi_bval', 'bids_pa_dwi',
                                       'bids_pa_dwi_bvec', 'bids_pa_dwi_bval'), list_tuple])

        parallelsubjects.inputs.bids_ap_dwi = bids_ap_dwi
        parallelsubjects.inputs.bids_ap_dwi_bvec = bids_ap_dwi_bvec
        parallelsubjects.inputs.bids_ap_dwi_bval = bids_ap_dwi_bval
        parallelsubjects.inputs.bids_pa_dwi = bids_pa_dwi
        parallelsubjects.inputs.bids_pa_dwi_bvec = bids_pa_dwi_bvec
        parallelsubjects.inputs.bids_pa_dwi_bval = bids_pa_dwi_bval

        self.connect([
            (parallelsubjects, self.input_node, [('bids_ap_dwi', 'bids_ap_dwi')]),
            (parallelsubjects, self.input_node, [('bids_ap_dwi_bvec', 'bids_ap_dwi_bvec')]),
            (parallelsubjects, self.input_node, [('bids_ap_dwi_bval', 'bids_ap_dwi_bval')]),
            (parallelsubjects, self.input_node, [('bids_pa_dwi', 'bids_pa_dwi')]),
            (parallelsubjects, self.input_node, [('bids_pa_dwi_bvec', 'bids_pa_dwi_bvec')]),
            (parallelsubjects, self.input_node, [('bids_pa_dwi_bval', 'bids_pa_dwi_bval')]),
            ])


    def build_output_node(self):
        """Build and connect an output node to the pipeline.
        """

        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        import nipype.interfaces.io as nio
        from os.path import join
        import dwi_preprocessing_noddi_utils as utils

        # Find container path from DWI filename
        # =====================================
        get_identifiers = npe.MapNode(nutil.Function(
            input_names=['in_file', 'caps_directory'], output_names=['base_directory', 'subst_tuple_list'],
            function=utils.get_subid_sesid), name='get_subid_sesid', iterfield=['in_file'])
        get_identifiers.inputs.caps_directory = self.caps_directory

        # datasink
        datasink = npe.MapNode(nio.DataSink(infields=['@preproc_bval', '@preproc_dwi', '@preproc_bvec', '@b0_mask']),
                              name='datasinker',
                              iterfield=['base_directory', 'substitutions', '@preproc_bval', '@preproc_dwi', '@preproc_bvec', '@b0_mask'])
        self.connect([
            # datasink
            (self.input_node, get_identifiers, [('bids_ap_dwi', 'in_file')]),
            (get_identifiers, datasink, [('base_directory', 'base_directory')]),
            (get_identifiers, datasink, [('subst_tuple_list', 'substitutions')]),
            # original files
            (self.output_node, datasink, [('preproc_bval', '@preproc_bval')]),
            (self.output_node, datasink, [('preproc_dwi', '@preproc_dwi')]),
            (self.output_node, datasink, [('preproc_bvec', '@preproc_bvec')]),
            (self.output_node, datasink, [('b0_mask', '@b0_mask')])
        ])

    def build_core_nodes(self):
        """Build and connect the core nodes of the pipeline.
        """
        import dwi_preprocessing_noddi_utils as utils

        one_subj_noddi = utils.noddi_preprocessing_twoped(self.caps_directory, epi_params=self.parameters['epi_param'],
                                                    alt_epi_params=self.parameters['alt_epi_params'])

        # Connection
        # ==========
        self.connect([
            (self.input_node,      one_subj_noddi,    [('bids_ap_dwi',    'inputnode.in_file')]),
            (self.input_node,      one_subj_noddi,    [('bids_ap_dwi_bvec',    'inputnode.in_bvec')]),
            (self.input_node,      one_subj_noddi,    [('bids_ap_dwi_bval',    'inputnode.in_bval')]),
            (self.input_node,      one_subj_noddi,    [('bids_pa_dwi',    'inputnode.alt_file')]),
            (self.input_node,      one_subj_noddi,    [('bids_pa_dwi_bvec',    'inputnode.alt_bvec')]),
            (self.input_node,      one_subj_noddi,    [('bids_pa_dwi_bval',    'inputnode.alt_bval')]),

            # output
            (one_subj_noddi, self.output_node, [('outputnode.ecc_out_file', 'preproc_dwi')]),  # noqa
            (one_subj_noddi, self.output_node, [('outputnode.out_bvec', 'preproc_bvec')]),  # noqa
            (one_subj_noddi, self.output_node, [('outputnode.original_merged_bval', 'preproc_bval')]),  # noqa
            (one_subj_noddi, self.output_node, [('outputnode.out_mask', 'b0_mask')])  # noqa
        ])
