"""dwi_processing_noddi - Clinica Pipeline.
This file has been generated automatically by the `clinica generate template`
command line tool. See here for more details: https://gitlab.icm-institute.org/aramislab/clinica/wikis/docs/InteractingWithClinica.
"""

# WARNING: Don't put any import statement here except if it's absolutly
# necessary. Put it *inside* the different methods.
# Otherwise it will slow down the dynamic loading of the pipelines list by the
# command line tool.
import clinica.pipelines.engine as cpe


class DwiProcessingNoddi(cpe.Pipeline):
    """dwi_processing_noddi SHORT DESCRIPTION.

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
        A clinica pipeline object containing the dwi_processing_noddi pipeline.

    Raises:


    Example:
        >>> import DwiProcessingNoddi
        >>> pipeline = dwi_processing_noddi('~/MYDATASET_BIDS', '~/MYDATASET_CAPS')
        >>> pipeline.parameters = {
        >>>     # ...
        >>> }
        >>> pipeline.base_dir = '/tmp/'
        >>> pipeline.run()
    """


    def check_custom_dependencies(self):
        """Check dependencies that can not be listed in the `info.json` file.
        """
        from clinica.utils.check_dependency import check_noddi_matlab_toolbox, check_nifti_matlib_toolbox
        _ = check_noddi_matlab_toolbox()
        _ = check_nifti_matlib_toolbox()

    def get_input_fields(self):
        """Specify the list of possible inputs of this pipeline.

        Returns:
            A list of (string) input fields name.
        """

        return ['subject_id_list', 
                'noddi_preprocessed_dwi', 
                'noddi_preprocessed_bvec', 
                'noddi_preprocessed_bval',
                'noddi_preprocessed_mask', 
                'n_procs', 
                'noddi_toolbox_dir', 
                'nifti_matlib_dir']  # Fill here the list


    def get_output_fields(self):
        """Specify the list of possible outputs of this pipeline.

        Returns:
            A list of (string) output fields name.
        """

        return ['fit_icvf', 'fit_isovf', 'fit_od']  # Fill here the list


    def build_input_node(self):
        """Build and connect an input node to the pipeline.
        """

        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        import dwi_processing_noddi_utils as utils

        subject_id_list, noddi_preprocessed_dwi, noddi_preprocessed_bvec, noddi_preprocessed_bval, noddi_preprocessed_mask = utils.grab_noddi_preprocessed_files(
            self.caps_directory, self.tsv_file)

        read_parameters_node = npe.Node(name="LoadingCLIArguments",
                                        interface=nutil.IdentityInterface(
                                            fields=self.get_input_fields(),
                                            mandatory_inputs=True))
        read_parameters_node.inputs.bvalue_str = self.parameters['n_procs']['n_procs']
        read_parameters_node.inputs.subject_id_list = subject_id_list
        read_parameters_node.inputs.noddi_preprocessed_dwi = noddi_preprocessed_dwi
        read_parameters_node.inputs.noddi_preprocessed_bvec = noddi_preprocessed_bvec
        read_parameters_node.inputs.noddi_preprocessed_bval = noddi_preprocessed_bval
        read_parameters_node.inputs.noddi_preprocessed_mask = noddi_preprocessed_mask
        read_parameters_node.inputs.noddi_toolbox_dir = self.parameters['noddi_toolbox_dir']['noddi_toolbox_dir']
        read_parameters_node.inputs.nifti_matlib_dir = self.parameters['nifti_matlib_dir']['nifti_matlib_dir']

        self.connect([
            (read_parameters_node,      self.input_node,    [('subject_id_list',    'subject_id_list')]),
            (read_parameters_node,      self.input_node,    [('noddi_preprocessed_dwi',    'noddi_preprocessed_dwi')]),
            (read_parameters_node,      self.input_node,    [('noddi_preprocessed_bvec',    'noddi_preprocessed_bvec')]),
            (read_parameters_node,      self.input_node,    [('noddi_preprocessed_bval',    'noddi_preprocessed_bval')]),
            (read_parameters_node,      self.input_node,    [('noddi_preprocessed_mask',    'noddi_preprocessed_mask')]),
            (read_parameters_node,      self.input_node,    [('n_procs',    'n_procs')]),
            (read_parameters_node,      self.input_node,    [('noddi_toolbox_dir',    'noddi_toolbox_dir')]),
            (read_parameters_node,      self.input_node,    [('nifti_matlib_dir',    'nifti_matlib_dir')]),
        ])


    def build_output_node(self):
        """Build and connect an output node to the pipeline.
        """
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        import nipype.interfaces.io as nio
        import dwi_processing_noddi_utils as utils

        # Find container path from DWI filename
        # =====================================
        get_identifiers = npe.MapNode(nutil.Function(
            input_names=['in_file', 'caps_directory'], output_names=['base_directory', 'subst_tuple_list'],
            function=utils.get_subid_sesid), name='get_subid_sesid', iterfield=['in_file'])
        get_identifiers.inputs.caps_directory = self.caps_directory

        ### datasink
        datasink = npe.MapNode(nio.DataSink(infields=['@fit_icvf', '@fit_isovf', '@fit_od']),
                              name='datasinker',
                              iterfield=['base_directory', 'substitutions', '@fit_icvf', '@fit_isovf', '@fit_od'])

        self.connect([
            ### datasink
            (self.output_node, get_identifiers, [('fit_icvf', 'in_file')]),
            (get_identifiers, datasink, [('base_directory', 'base_directory')]),
            (get_identifiers, datasink, [('subst_tuple_list', 'substitutions')]),
            ## original files
            (self.output_node, datasink, [('fit_icvf', '@fit_icvf')]),
            (self.output_node, datasink, [('fit_isovf', '@fit_isovf')]),
            (self.output_node, datasink, [('fit_od', '@fit_od')]),
        ])

    def build_core_nodes(self):
        """Build and connect the core nodes of the pipeline.
        """

        import dwi_processing_noddi_utils as utils
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe

        processing_pipeline = utils.matlab_noddi_processing(self.caps_directory,
                                                            num_cores=self.parameters['n_procs']['n_procs'], bStep=self.parameters['bvalue_str']['bvalue_str'])

        # Connection
        # ==========
        self.connect([
            (self.input_node,      processing_pipeline,    [('subject_id_list',    'inputnode.subject_id_list')]),
            (self.input_node,      processing_pipeline,    [('noddi_preprocessed_dwi',    'inputnode.noddi_preprocessed_dwi')]),
            (self.input_node,      processing_pipeline,    [('noddi_preprocessed_bvec',    'inputnode.noddi_preprocessed_bvec')]),
            (self.input_node,      processing_pipeline,    [('noddi_preprocessed_bval',    'inputnode.noddi_preprocessed_bval')]),
            (self.input_node,      processing_pipeline,    [('noddi_preprocessed_mask',    'inputnode.noddi_preprocessed_mask')]),
            (self.input_node,      processing_pipeline,    [('noddi_toolbox_dir',    'inputnode.noddi_toolbox_dir')]),
            (self.input_node,      processing_pipeline,    [('nifti_matlib_dir',    'inputnode.nifti_matlib_dir')]),
            ## output
            (processing_pipeline, self.output_node, [('outputnode.fit_icvf', 'fit_icvf')]),  # noqa
            (processing_pipeline, self.output_node, [('outputnode.fit_isovf', 'fit_isovf')]),  # noqa
            (processing_pipeline, self.output_node, [('outputnode.fit_od', 'fit_od')])  # noqa
        ])
