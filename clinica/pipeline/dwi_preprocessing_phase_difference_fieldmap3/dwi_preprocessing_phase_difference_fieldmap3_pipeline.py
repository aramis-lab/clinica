"""Dwi Preprocessing Phase Difference Fieldmap3 - Clinica Pipeline.
This file has been generated automatically by the `clinica generate template`
command line tool. See here for more details: https://gitlab.icm-institute.org/aramis/clinica/wikis/docs/InteractingWithClinica.
"""

# WARNING: Don't put any import statement here except if it's absolutly
# necessary. Put it *inside* the different methods.
# Otherwise it will slow down the dynamic loading of the pipelines list by the
# command line tool.
import clinica.pipeline.engine as cpe


class DwiPreprocessingPhaseDifferenceFieldmap3(cpe.Pipeline):
    """Dwi Preprocessing Phase Difference Fieldmap3 SHORT DESCRIPTION.

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
        A clinica pipeline object containing the Dwi Preprocessing Phase Difference Fieldmap3 pipeline.

    Raises:


    Example:
        >>> from dwi_preprocessing_phase_difference_fieldmap3 import DwiPreprocessingPhaseDifferenceFieldmap3
        >>> pipeline = FMRIPreprocessing('~/MYDATASET_BIDS', '~/MYDATASET_CAPS')
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

        return ['in_dwi', 'in_fmap_magnitude', 'in_fmap_phasediff', 'in_bval', 'in_bvec', 'work_dir', 'register_fmap_on_b0', 'register_b0_on_b0', 'echospacing', 'delta_te', 'enc_dir']


    def get_output_fields(self):
        """Specify the list of possible outputs of this pipeline.

        Returns:
            A list of (string) output fields name.
        """

        return ['out_bvecs', 'out_bvals', 'out_preprocessed_dwi', 'out_b0_mask']


    def build_input_node(self):
        """Build and connect an input node to the pipeline.
        """


        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe


        # This node is supposedly used to load BIDS inputs when this pipeline is
        # not already connected to the output of a previous Clinica pipeline.
        # For the purpose of the example, we simply read input arguments given
        # by the command line interface and transmitted here through the
        # `self.parameters` dictionary and pass it to the `self.input_node` to
        # further by used as input of the core nodes.

        read_parameters_node = npe.Node(name="LoadingCLIArguments",
                                        interface=nutil.IdentityInterface(
                                            fields=self.get_input_fields(),
                                            mandatory_inputs=True), synchronize=True)

        read_parameters_node.inputs.register_fmap_on_b0 = self.parameters['register_fmap_on_b0']
        read_parameters_node.inputs.register_b0_on_b0 = self.parameters['register_b0_on_b0']
        read_parameters_node.inputs.delta_te = self.parameters['delta_te']
        read_parameters_node.inputs.echospacing = self.parameters['echospacing']
        read_parameters_node.inputs.enc_dir = self.parameters['enc_dir']
        read_parameters_node.inputs.in_dwi = self.bids_layout.get(return_type='file', type='dwi', extensions='nii|nii.gz', session='|'.join(self.sessions).replace('ses-',''),
                             subject='|'.join(self.subjects).replace('sub-',''))
        read_parameters_node.inputs.in_fmap_magnitude = self.bids_layout.get(return_type='file', type='magnitude1', extensions='nii|nii.gz', session='|'.join(self.sessions).replace('ses-',''),
                             subject='|'.join(self.subjects).replace('sub-',''))
        read_parameters_node.inputs.in_fmap_phasediff = self.bids_layout.get(return_type='file', type='phasediff', extensions='nii|nii.gz', session='|'.join(self.sessions).replace('ses-',''),
                             subject='|'.join(self.subjects).replace('sub-',''))
        read_parameters_node.inputs.in_bval = self.bids_layout.get(return_type='file', type='dwi', extensions='bval', session='|'.join(self.sessions).replace('ses-',''),
                             subject='|'.join(self.subjects).replace('sub-',''))
        read_parameters_node.inputs.in_bvec = self.bids_layout.get(return_type='file', type='dwi', extensions='bvec', session='|'.join(self.sessions).replace('ses-',''),
                             subject='|'.join(self.subjects).replace('sub-',''))
        read_parameters_node.inputs.work_dir = [self.base_dir + '/DwiPreprocessingPhaseDifferenceFieldmap3/' + self.subjects[i] for i in xrange(len(self.subjects))]

        list_tuple = []
        for i in xrange(len(self.subjects)):
            ls = []
            ls.append(read_parameters_node.inputs.in_dwi[i])
            ls.append(read_parameters_node.inputs.in_bval[i])
            ls.append(read_parameters_node.inputs.in_bvec[i])
            ls.append(read_parameters_node.inputs.in_fmap_magnitude[i])
            ls.append(read_parameters_node.inputs.in_fmap_phasediff[i])
            ls.append(read_parameters_node.inputs.work_dir[i])
            ls = tuple(ls)
            list_tuple.append(ls)

        read_parameters_node.iterables = [('in_dwi', 'in_bval', 'in_bvec', 'in_fmap_magnitude', 'in_fmap_phasediff', 'work_dir'),
                                          list_tuple]

        self.connect([(read_parameters_node,      self.input_node,    [('in_dwi',    'in_dwi')])])
        self.connect([(read_parameters_node,      self.input_node,    [('in_fmap_magnitude',    'in_fmap_magnitude')])])
        self.connect([(read_parameters_node,      self.input_node,    [('in_fmap_phasediff',    'in_fmap_phasediff')])])
        self.connect([(read_parameters_node,      self.input_node,    [('in_bval',    'in_bval')])])
        self.connect([(read_parameters_node,      self.input_node,    [('in_bvec',    'in_bvec')])])
        self.connect([(read_parameters_node,      self.input_node,    [('register_fmap_on_b0',    'register_fmap_on_b0')])])
        self.connect([(read_parameters_node,      self.input_node,    [('register_b0_on_b0',    'register_b0_on_b0')])])
        self.connect([(read_parameters_node,      self.input_node,    [('echospacing',    'echospacing')])])
        self.connect([(read_parameters_node,      self.input_node,    [('delta_te',    'delta_te')])])
        self.connect([(read_parameters_node,      self.input_node,    [('enc_dir',    'enc_dir')])])
        self.connect([(read_parameters_node,      self.input_node,    [('work_dir',    'work_dir')])])

    def build_output_node(self):
        """Build and connect an output node to the pipeline.
        """

        # In the same idea as the input node, this output node is supposedly
        # used to write the output fields in a CAPS. It should be executed only
        # if this pipeline output is not already connected to a next Clinica
        # pipeline.

        pass


    def build_core_nodes(self):
        """Build and connect the core nodes of the pipeline.
        """


        import dwi_preprocessing_phase_difference_fieldmap3_utils as utils
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        from clinica.utils.check_dependency import check_ants, check_fsl

        check_ants()
        check_fsl()

        full_pipe = npe.Node(nutil.Function(input_names=['in_dwi', 'in_bvals', 'in_bvecs', 'in_fmap_magnitude', 'in_fmap_phasediff', 'delta_te', 'echospacing', 'enc_dir',
                'register_fmap_on_b0', 'register_b0_on_b0', 'work_dir', 'n_procs', 'caps_dir'],
                                           output_names=[],
                                           function=utils.diffusion_preprocessing_phasediff_fieldmap),
                              name='full_pipe_preprocessing'

        )
        full_pipe.inputs.register_fmap_on_b0 = self.parameters['register_fmap_on_b0']
        full_pipe.inputs.register_b0_on_b0 = self.parameters['register_b0_on_b0']
        full_pipe.inputs.work_dir = self.base_dir
        full_pipe.inputs.delta_te = self.parameters['delta_te']
        full_pipe.inputs.echospacing = self.parameters['echospacing']
        full_pipe.inputs.enc_dir = self.parameters['enc_dir']
        full_pipe.inputs.n_procs = self.parameters['n_procs']
        full_pipe.inputs.caps_dir = self.caps_directory

        self.connect([
            (self.input_node, full_pipe, [('in_dwi', 'in_dwi', ),
                              ('in_bval', 'in_bvals'),
                              ('in_bvec', 'in_bvecs')]),
            (self.input_node, full_pipe, [('in_fmap_magnitude', 'in_fmap_magnitude')]),
            (self.input_node, full_pipe, [('in_fmap_phasediff', 'in_fmap_phasediff')]),
            (self.input_node, full_pipe, [('work_dir', 'work_dir')]),

            # Outputnode:
            # (full_pipe, self.output_node, [('outputnode.out_bvecs', 'out_bvecs')]),
            # (full_pipe, self.output_node, [('outputnode.out_bvals', 'out_bvals')]),
            # (full_pipe, self.output_node, [('outputnode.out_file', 'out_preprocessed_dwi')]),
            # (full_pipe, self.output_node, [('outputnode.b0_mask', 'out_b0_mask')])

        ])
        return self